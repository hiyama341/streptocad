import os
import pytest
from src.utils import list_of_objects_in_a_dir

# Test the normal behavior with a temporary directory
def test_list_of_objects_in_a_dir_with_files(tmp_path):
    # Create files and directories in the temporary path
    (tmp_path / "file1.txt").write_text("content")
    (tmp_path / "file2.txt").write_text("content")
    subdir = tmp_path / "subdir"
    subdir.mkdir()  # Corrected this line to create a directory without mode as a string
    subdir.joinpath("file3.txt").write_text("content")

    # Call the function with the path to the temp directory
    files_list = list_of_objects_in_a_dir(tmp_path)

    # Check if the function returns only files, not directories
    assert len(files_list) == 2
    assert "file1.txt" in files_list
    assert "file2.txt" in files_list
    assert "file3.txt" not in files_list  # file3.txt is in a subdir


# Test the function with an empty directory
def test_list_of_objects_in_a_dir_empty(tmp_path):
    files_list = list_of_objects_in_a_dir(tmp_path)
    assert files_list == []  # Should return an empty list

# Test the function with a non-existent directory
def test_list_of_objects_in_a_dir_nonexistent(tmp_path):
    # Define a path that does not exist
    non_existent_path = tmp_path / "nonexistent_dir"
    
    # We expect a FileNotFoundError
    with pytest.raises(FileNotFoundError):
        list_of_objects_in_a_dir(non_existent_path)

