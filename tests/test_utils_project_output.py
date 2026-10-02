"""Tests for the output-building helpers in streptocad.utils.

These assemble the project package every workflow hands the user, so they cover
each content type the writers branch on: lists of Dseqrecords, single records,
DataFrames and plain strings.
"""

import io
import json
import os
import zipfile

import pandas as pd
import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pydna.dseqrecord import Dseqrecord

from streptocad.utils import (
    ProjectDirectory,
    create_primer_df_from_dict,
    dataframe_to_seqrecords,
    extract_metadata_to_dataframe,
    generate_header,
    generate_project_directory_structure,
    get_top_sgRNAs,
)
from tests.conftest import stub_html_exporter

# --- dataframe_to_seqrecords -------------------------------------------


def test_dataframe_to_seqrecords_labels_each_sgrna():
    df = pd.DataFrame(
        {
            "locus_tag": ["SCO1", "SCO2"],
            "sgrna_loc": [10, 20],
            "sgrna": ["ATCGATCGATCGATCGATCG", "CGTACGTACGTACGTACGTA"],
        }
    )

    records = dataframe_to_seqrecords(df)

    assert [r.id for r in records] == ["SCO1_10", "SCO2_20"]
    assert str(records[0].seq) == "ATCGATCGATCGATCGATCG"
    for record in records:
        (feature,) = record.features
        assert feature.type == "sgRNA"
        assert feature.qualifiers["label"] == record.id
        assert len(feature.location) == len(record.seq)


def test_dataframe_to_seqrecords_on_an_empty_frame():
    df = pd.DataFrame(columns=["locus_tag", "sgrna_loc", "sgrna"])
    assert dataframe_to_seqrecords(df) == []


# --- get_top_sgRNAs ----------------------------------------------------


def test_get_top_sgrnas_keeps_the_fewest_off_targets_per_locus():
    df = pd.DataFrame(
        {
            "locus_tag": ["A", "A", "A", "B", "B"],
            "sgrna": ["a1", "a2", "a3", "b1", "b2"],
            "off_target_n": [5, 1, 3, 9, 2],
        }
    )

    top = get_top_sgRNAs(df, number_of_sgRNAs=2)

    by_locus = top.groupby("locus_tag")["sgrna"].apply(set).to_dict()
    assert by_locus == {"A": {"a2", "a3"}, "B": {"b1", "b2"}}


def test_get_top_sgrnas_asking_for_more_than_exist_returns_all():
    df = pd.DataFrame(
        {"locus_tag": ["A", "A"], "sgrna": ["a1", "a2"], "off_target_n": [1, 2]}
    )
    assert len(get_top_sgRNAs(df, number_of_sgRNAs=10)) == 2


# --- create_primer_df_from_dict ----------------------------------------


def primer_record(gene):
    return {
        "gene_name": gene,
        "up_fwd_p_anneal": Seq("ATCG"),
        "up_reverse_p_anneal": Seq("CGTA"),
        "tm_up_forwar_p": 60.0,
        "tm_up_reverse_p": 61.0,
        "ta_up": 55.0,
        "up_forwar_primer_str": "ATCGATCG",
        "up_reverse_primer_str": "CGTACGTA",
        "up_forwar_p_name": f"{gene}_up_fwd",
        "up_reverse_p_name": f"{gene}_up_rev",
        "dw_fwd_p_anneal": Seq("GCTA"),
        "dw_reverse_p_anneal": Seq("TAGC"),
        "tm_dw_forwar_p": 58.0,
        "tm_dw_reverse_p": 59.0,
        "ta_dw": 54.0,
        "dw_forwar_primer_str": "GCTAGCTA",
        "dw_reverse_primer_str": "TAGCTAGC",
        "dw_forwar_p_name": f"{gene}_dw_fwd",
        "dw_reverse_p_name": f"{gene}_dw_rev",
    }


def test_create_primer_df_makes_two_rows_per_gene():
    df = create_primer_df_from_dict([primer_record("geneA"), primer_record("geneB")])

    assert len(df) == 4
    assert list(df["direction"]) == ["upstream", "downstream"] * 2
    assert set(df["template"]) == {"geneA", "geneB"}

    upstream = df[df["direction"] == "upstream"].iloc[0]
    assert upstream["f_primer_anneal(5-3)"] == "ATCG"
    assert upstream["f_primer_anneal_length"] == 4
    assert upstream["f_primer_length"] == len("ATCGATCG")
    assert upstream["ta"] == 55.0
    assert upstream["f_primer_name"] == "geneA_up_fwd"


def test_create_primer_df_on_no_records_keeps_its_columns():
    df = create_primer_df_from_dict([])
    assert df.empty
    assert "f_primer_anneal(5-3)" in df.columns


# --- extract_metadata_to_dataframe -------------------------------------


@pytest.fixture
def assembled():
    plasmid = Dseqrecord(Seq("ATGC" * 30), circular=True)
    plasmid.name = "pBase"
    plasmid.id = "pBase"
    built = []
    for i in range(2):
        record = Dseqrecord(Seq("ATGC" * 40), circular=True)
        record.name = f"pNew_{i}"
        built.append(record)
    return plasmid, built


def test_extract_metadata_describes_each_plasmid(assembled):
    plasmid, built = assembled

    df = extract_metadata_to_dataframe(built, plasmid, ["sgRNA_1", "sgRNA_2"])

    assert len(df) == 2
    assert list(df["integration"]) == ["sgRNA_1", "sgRNA_2"]
    assert set(df["original_plasmid"]) == {"pBase"}
    assert list(df["size"]) == [len(r.seq) for r in built]
    assert df["date"].notna().all()


# --- generate_header ---------------------------------------------------


def test_generate_header_counts_and_appends_the_capture():
    header = generate_header(["p1", "p2"], ["s1"], [1, 2, 3], "captured text\n")

    assert "generated 2 plasmids from 1 sequences" in header
    assert "generated 3 primers" in header
    assert header.endswith("captured text\n")


def test_generate_header_with_no_capture():
    assert generate_header([], [], [], "").endswith("\n\n\n")


# --- ProjectDirectory --------------------------------------------------


@pytest.fixture
def package_inputs():
    record = Dseqrecord(Seq("ATGCATGCATGC"), circular=True)
    record.name = "single"
    record.id = "single"

    seq_list = []
    for i in range(2):
        item = Dseqrecord(Seq("GGCCGGCCGGCC"))
        item.name = f"listed_{i}"
        item.id = f"listed_{i}"
        seq_list.append(item)

    return {
        "input_files": [
            {"name": "many.gb", "content": seq_list},
            {"name": "one.gb", "content": record},
        ],
        "output_files": [
            {"name": "table.csv", "content": pd.DataFrame({"a": [1, 2]})},
            {"name": "notes.log", "content": "a log line\n"},
            {"name": "built.gb", "content": seq_list},
        ],
        "input_values": {"settings": {"tm": 65}},
    }


def read_zip(archive):
    return zipfile.ZipFile(io.BytesIO(archive))


def test_project_directory_packages_every_content_type(package_inputs):
    project = ProjectDirectory(project_name="demo", **package_inputs)

    with stub_html_exporter():
        archive = project.create_directory_structure()

    with read_zip(archive) as zf:
        names = zf.namelist()
        basenames = {n.rsplit("/", 1)[-1] for n in names}

        # On the input side a list of records becomes one multi-record file...
        assert {"many.gb", "one.gb", "input_values.json"} <= basenames
        # ...while on the output side it is split, named after each record.
        assert {"listed_0_0.gb", "listed_1_1.gb"} <= basenames
        assert {"table.csv", "notes.log"} <= basenames

        assert (
            zf.read(next(n for n in names if n.endswith("notes.log")))
            == b"a log line\n"
        )
        assert (
            "a,b"
            not in zf.read(next(n for n in names if n.endswith("table.csv"))).decode()
        )

        values = json.loads(zf.read(next(n for n in names if n.endswith(".json"))))
        assert values == {"settings": {"tm": 65}}

        assert all(n.startswith("demo/") for n in names)
        assert any(n == "demo/inputs/" for n in names)
        assert any(n == "demo/outputs/" for n in names)


def test_project_directory_can_skip_the_directory_entries(package_inputs):
    project = ProjectDirectory(project_name="demo", **package_inputs)

    with stub_html_exporter():
        archive = project.create_directory_structure(create_directories=False)

    with read_zip(archive) as zf:
        assert "demo/inputs/" not in zf.namelist()
        assert any(n.endswith("one.gb") for n in zf.namelist())


def test_project_directory_renders_markdown_protocols(package_inputs, tmp_path):
    protocol = tmp_path / "protocol.md"
    protocol.write_text("# Step one\n")

    project = ProjectDirectory(
        project_name="demo", markdown_file_paths=[str(protocol)], **package_inputs
    )
    with stub_html_exporter():
        archive = project.create_directory_structure()

    with read_zip(archive) as zf:
        html_name = next(n for n in zf.namelist() if n.endswith("protocol.html"))
        assert "Step one" in zf.read(html_name).decode()
        assert not any(n.endswith("protocol.md") for n in zf.namelist())


def test_project_directory_records_what_it_wrote(package_inputs):
    project = ProjectDirectory(project_name="demo", **package_inputs)
    with stub_html_exporter():
        project.create_directory_structure()

    assert project.project_dir_structure
    assert any(path.endswith("table.csv") for path in project.project_dir_structure)


def test_get_zip_file_rewinds_the_buffer(package_inputs):
    project = ProjectDirectory(project_name="demo", **package_inputs)
    with stub_html_exporter():
        archive = project.create_directory_structure()

    buffer = project.get_zip_file()
    assert buffer.tell() == 0
    assert buffer.read() == archive


# --- generate_project_directory_structure ------------------------------


def test_generate_project_directory_structure_writes_to_disk(package_inputs, tmp_path):
    cwd = os.getcwd()
    os.chdir(tmp_path)
    try:
        structure = generate_project_directory_structure(
            project_name="on_disk", **package_inputs
        )
    finally:
        os.chdir(cwd)

    assert structure
    written = {p for p in structure if (tmp_path / p).exists()}
    assert written, "nothing was written to disk"
    assert (tmp_path / "on_disk" / "inputs").is_dir()
    assert (tmp_path / "on_disk" / "outputs").is_dir()


def test_generate_project_directory_structure_can_plan_without_writing(
    package_inputs, tmp_path
):
    cwd = os.getcwd()
    os.chdir(tmp_path)
    try:
        structure = generate_project_directory_structure(
            project_name="planned", create_directories=False, **package_inputs
        )
    finally:
        os.chdir(cwd)

    assert structure
    assert not (tmp_path / "planned").exists()
