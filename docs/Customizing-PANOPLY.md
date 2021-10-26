### Customizing PANOPLY

PANOPLY can be customized using the flexibility and customization abilities of [Terra](https://app.terra.bio/). [Terra Support] (https://support.terra.bio/) includes documentation, tutorials and guides to help users add new tasks and customize workflows. Customizing PANOPLY will require some programming ability and comfort with command-line environments. Here, we present links to Terra support/documentation to point users to the appropriate resources.

#### Customizing Workflows
PANOPLY workflows can be customized by (i) adding/removing tasks from existing workflows or (ii) creating new workflows combining existing tasks in different ways. The following provide pointers to how tasks and workflows are structured. 

* [Introduction to WDL](https://support.terra.bio/hc/en-us/articles/360037117492-Getting-Started-with-WDL)
* [Create, edit, and share a new workflow](https://support.terra.bio/hc/en-us/articles/360031366091-Create-edit-and-share-a-new-workflow)
* [Multi-step wrokflows](https://support.terra.bio/hc/en-us/articles/360037128072--2-howto-Write-a-simple-multi-step-workflow)

#### Adding new tasks
Users can add new algorithms to PANOPLY by encapsulating them into a `docker`, and then using these dockers in tasks and workflows. 

* [Create a docker image](https://support.terra.bio/hc/en-us/articles/360024737591-Make-a-Docker-container-image-the-easy-way-using-a-base-image)

Users can explore the `hydrant` and `src` subdirectories in the [PANOPLY repository](https://github.com/broadinstitute/PANOPLY) to see examples of how docker images are created for R code, how these images are converted to tasks using WDL, and how multiple tasks are strung together into a full pipeline using WDL workflows.