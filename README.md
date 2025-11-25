
# Nextflow 101 Training Repository

## Purpose

This repository provides a collection of step-by-step Nextflow scripts used during the **Nextflow 101 training session**. It is designed to support learners with a software development background (e.g., C++, JavaScript) who are new to Nextflow. The training covers fundamental concepts, best practices, and hands-on pipeline creation.

By following the script versions in this repo, trainees can:

* Learn the basics of Nextflow DSL2
* Understand input/output channels and process definitions
* Run and modify example pipelines locally or on a cloud platform
* Progress at their own pace by either writing code from scratch (starting at the `start` tag) or jumping to a later tag to see a completed solution

## Using Git Tags for Navigation

To help learners follow along during and after the training, **Git tags** have been used to mark important milestones in the development of the pipeline.

### What is a Git tag?

A Git tag is a label that points to a specific commit in the Git history. Tags are often used to mark release points or important stages of a project.

In this training repo, we use tags to mark key milestones in pipeline development:

* `start` — the initial scaffold with minimal Nextflow configuration (recommended starting point for writing code yourself)
* `fastqc_only` — pipeline with FastQC analysis implemented
* `fastp_and_fastqc` — pipeline with both FastQC and fastp steps
* `multiqc` — final version with MultiQC report generation

### How to list available tags

```bash
git tag
```

### How to checkout a specific tag (in detached HEAD mode)

```bash
git checkout tags/<tag_name>
```

For example, to jump to the FastQC-only solution:

```bash
git checkout tags/fastqc_only
```

Alternatively, to view the code at a tag without leaving your branch:

```bash
git show <tag_name>
```

### How to compare versions using `git diff`

To see what changed between two versions:

```bash
git diff start fastqc_only
git diff fastqc_only fastp_and_fastqc
git diff fastp_and_fastqc multiqc
```

This will highlight what was added, removed, or modified between the stages.

## Additional Tips for Beginners

* **Stay on a branch**: After exploring a tag, return to your working branch with:

  ```bash
  git checkout main
  ```

* **Clone the repo with history**: Make sure to clone with full history (not shallow) so tags are available:

  ```bash
  git clone <repo-url>
  cd nextflow-101
  ```

* **Use Visual Tools**: If you're using VS Code or GitHub Desktop, tags can be accessed from the GUI under the version history panel.

* **Tag your own work**: As you experiment, use tags to bookmark your progress:

  ```bash
  git tag my_first_pipeline
  git push origin my_first_pipeline
  ```

## Getting Started

1. Install [Nextflow](https://www.nextflow.io/)
2. Checkout the `start` tag and follow the guided steps in `main.nf` to write the code yourself, **or**
3. Jump to any later tag (`fastqc_only`, `fastp_and_fastqc`, `multiqc`) to see a completed solution for that stage
4. Use `git diff` to compare your work with the provided solutions

We hope this structure helps you learn Nextflow more effectively. Happy scripting!
