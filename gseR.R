library(googleComputeEngineR)
project = "scmerge"
# zone = "us-east1-b"
zone = "australia-southeast1-a"

gce_global_project(project)
gce_global_zone(zone)
# gce_get_project()
# gce_list_zones(project)
# gce_list_machinetype()$items

(tag = gce_tag_container("singlecellplus_docker"))

vm <- gce_vm(template = "rstudio", 
             name = "scp-name", 
             predefined_type = "n1-highmem-16",
             dynamic_image = tag,
             user = "rstudio", 
             password = "pushu")

# bash: gcloud compute ssh scp-name
vm <- gce_ssh_setup(vm,
                    username = "rstudio", 
                    key.pub = "~/.ssh/id_rsa.pub",
                    key.private = "~/.ssh/id_rsa")
# gce_ssh(vm, "echo foo", username = "rstudio")
gce_rstudio_adduser(instance = vm, username = "jean", password = "jean")
gce_rstudio_adduser(instance = vm, username = "kevin", password = "kevin")
gce_rstudio_adduser(instance = vm, username = "pushu", password = "pushu")

gce_vm_stop(instances = vm)
