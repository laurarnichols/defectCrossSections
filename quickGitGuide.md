# Quick Git Guide
## Connecting to Git Using an SSH Key
* Log in to the machine you want the local copy of the repo to live on
* Generate an SSH key using `ssh-keygen`, hitting enter at each of the prompts
* Your public SSH key will be stored in `~/.ssh/id_rsa.pub`
* Log in to your profile on a web browser and go to settings
* Select "SSH and GPG keys"
* Click new SSH key and copy and paste the contents of the `~/.ssh/id_rsa.pub` file into the large text box, adding a name to the smaller box above, if desired
* You can now clone a git repo from this machine using ssh (see [Forking and Cloning the Repo and Compiling](#forking-and-cloning-the-repo-and-compiling))

## Forking and Cloning Repo
* First, go the the main repo on a web browser and fork it using the fork button at the top right
* You will now have a copy of the repo under your username
* Clone the repo on the machine of your choice using ssh (must have ssh key set up; see [Connecting to Git Using an SSH Key](#connecting-to-git-using-an-ssh-key)):
```
git clone git@github.com:<your git username>/defectCrossSections.git
```
* You now have a local copy of your fork on your machine
