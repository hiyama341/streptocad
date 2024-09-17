# StreptoCAD

<img src="web_app/assets/StreptoCAD_logo Medium.jpeg" alt="StreptoCAD" width="200">

**StreptoCAD** is an open-source software toolbox designed to automate and streamline genome engineering in Streptomyces. This tool supports various CRISPR-based techniques and gene overexpression methods, significantly simplifying the genetic engineering process.

## Features

- **Automated Primer and sgRNA Design:** Automatically generates necessary DNA primers and sgRNA sequences for your target genes.
- **Plasmid Assembly Simulation:** Simulates plasmid assemblies and the resulting genomic modifications.
- **Six Design Workflows:** Supports workflows including overexpression library construction, base-editing, and in-frame deletions using CRISPR-Cas9 and CRISPR-Cas3 systems.
- **FAIR Compliance:** Ensures data is Findable, Accessible, Interoperable, and Reusable, promoting reproducibility and ease of data management.
- **User-Friendly:** Suitable for both experienced users and beginners, facilitating collaboration and standardized workflows.

## Why StreptoCAD?

Streptomyces is a prolific source of novel bioactive molecules, but current genetic engineering methods are inefficient and time-consuming. StreptoCAD addresses these challenges by automating the design process, reducing errors, and speeding up the development of genetically modified strains. This tool transforms complex genetic engineering tasks into straightforward, reproducible processes, enabling faster scientific advancements and discovery of new natural products.

## Workflows

StreptoCAD offers six distinct workflows for various genetic engineering tasks:

1. **[Overexpression Plasmid Library Construction](#workflow-1-overexpression-plasmid-library-construction):**

   - Can be used to overexpress target proteins - we experimentally validated this by overexprssing regulators.

2. **[Single CRISPR-BEST Plasmid Generation](#workflow-2-single-crispr-best-plasmid-generation):**

   - Base editing system in the genome of Streptomyces using single sgRNA for targeting.

3. **[Multiplexed CRISPR-BEST Plasmid Generation](#workflow-3-multiplexed-crispr-best-plasmid-generation):**

   - Multiplexed base-editing in the genome for high-throughput genetic studies.

4. **[CRISPRi Plasmid Generation](#workflow-4-crispri-plasmid-generation):**

   - Uses transcriptional interference to reversibly inactivate genes for functional studies.

5. **[CRISPR-Cas9](#workflow-5-in-frame-deletion-with-crispr-cas9):**

   - Can be used for random-sized or in-frame deletions with Cas9

6. **[CRISPR-Cas3](#workflow-6-in-frame-deletion-with-crispr-cas3):**
   - Can be used for random-sized or in-frame deletions with Cas3

## Experimental Validation

StreptoCAD's efficiency and user-friendliness were validated by designing and constructing overexpression strains in Streptomyces Göe40/10 in just eight weeks. This highlights the tool's capability to accelerate genome engineering projects.

## Future Developments

Future expansions will include additional genome engineering tools and integration with laboratory robotics systems for end-to-end automation, further enhancing the capabilities and efficiency of StreptoCAD.

## Get Started

Visit [www.streptocad.bioengineering.dtu.com](www.streptocad.bioengineering.dtu.com) to download StreptoCAD, access detailed documentation, and join the community of users and contributors working to advance Streptomyces research.

## Want to run StreptoCAD locally?

#### 1. Set up a Conda virtual environment (Why it's smart)

Using a Conda virtual environment is a great way to manage dependencies for your project. Conda makes it easy to create and manage isolated environments, ensuring your project’s libraries are kept separate from other projects and system-wide dependencies. This helps avoid compatibility issues and makes it simpler to reproduce your development environment.

To create a new Conda environment, run:

```bash
conda create --name streptocad python=3.11
```

Replace streptocad with your preferred environment name, and replace 3.11 with the specific version of Python you need.

Then activate it:

```bash
conda activate streptocad
```

#### 2. Install the requirements

Once your Conda environment is active, you can install the required dependencies from requirements.txt. This ensures your environment has all the necessary packages for the project. Use the following command:

```bash
pip install -r requirements.txt
```

(Note: Even though you're using Conda, pip is still used to install from requirements.txt.)

#### 3. Run the application

Finally, to run the StreptoCAD application, execute the following command:

```bash
python3 application.py
```

This will launch the application locally, and you're ready to go! Follow the url that your terminal shows.

## License

StreptoCAD is open-source and licensed under the MIT License.

## Contact

For questions or contributions, please contact [luclev@dtu.dk](mailto:luclev@dtu.dk).
