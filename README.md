# Nextflow 101 Training Repository

## Purpose

This repository provides a collection of step-by-step Nextflow scripts used during the **Nextflow 101 training session**. It is designed to support learners with a software development background (e.g., C++, JavaScript) who are new to Nextflow. The training covers fundamental concepts, best practices, and hands-on pipeline creation.

By following the script versions in this repo, trainees can:

* Learn the basics of Nextflow DSL2
* Understand input/output channels and process definitions
* Run and modify example pipelines locally or on a cloud platform

## Using Git Tags for Navigation

To help learners follow along during and after the training, **Git tags** have been used to mark important milestones in the development of the pipeline.

### What is a Git tag?

A Git tag is a label that points to a specific commit in the Git history. Tags are often used to mark release points or important stages of a project.

In this training repo, we use tags to mark:

* `start` — the initial scaffold with minimal Nextflow configuration
* `local_demo_done` — pipeline ready to be run locally
* `cloud_demo_done` — final version ready for cloud execution

### How to list available tags

```bash
git tag
```

### How to checkout a specific tag (in detached HEAD mode)

```bash
git checkout tags/local_demo_done
```

Alternatively, to explore without leaving your current branch:

```bash
git show local_demo_done
```

### How to compare versions using `git diff`

To see what changed between two versions:

```bash
git diff start local_demo_done
```

This will highlight what was added, removed, or modified between the two stages.

## Additional Tips for Beginners

* **Stay on a branch**: After exploring a tag, return to your working branch with:

  ```bash
  git checkout main
  ```

* **Clone the repo with history**: Make sure to clone with full history (not shallow) so tags are available:

  ```bash
  git clone https://github.com/your-org/nextflow-101.git
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
2. Checkout the `start` tag and follow the guided steps in `main.nf`
3. Progress through tags as the training evolves

We hope this structure helps you learn Nextflow more effectively. Happy scripting!
