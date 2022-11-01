# Quick Git Guide
## Connecting to Git Using an SSH Key
* Log in to the machine you want the local copy of the repo to live on
* Generate an SSH key using `ssh-keygen`, hitting enter at each of the prompts
* Your public SSH key will be stored in `~/.ssh/id_rsa.pub`
* Log in to your profile on a web browser and go to settings
* Select "SSH and GPG keys"
* Click new SSH key and copy and paste the contents of the `~/.ssh/id_rsa.pub` file into the large text box, adding a name to the smaller box above, if desired
* You can now clone a git repo from this machine using ssh

## Forking and Cloning Repo
* First, go the the main repo on a web browser and fork it using the fork button at the top right
* You will now have a copy of the repo under your username
* Clone the repo on the machine of your choice using ssh (must have ssh key set up; see [Connecting to Git Using an SSH Key](#connecting-to-git-using-an-ssh-key)):
```
git clone git@github.com:<your git username>/defectCrossSections.git
```
* You now have a local copy of your fork on your machine

## Contributing
1. Before you start a project, get an up-to-date version of the code (see [this Git documentation on merging from upstream repos](https://docs.github.com/en/github/collaborating-with-issues-and-pull-requests/merging-an-upstream-repository-into-your-fork))
2. To make changes, first create a new branch using `git checkout -b NAME-OF-BRANCH`. You can also switch to an existing branch
   using `git checkout NAME-OF-BRANCH`. Make the branch name represent the goal of your project (e.g., `parallelize-export`).
3. You can now make changes to files and git will track the changes. At any point, you can use `git status` to see what files you have changed.
4. A single project will have many tasks representing isolated changes. Each time you make an isolated change, create a commit:
	* Stage the file or folder using `git add FILE OR FOLDER`. 
	* Complete commit using `git commit -m "COMMENT TO DESCRIBE THE INTENTION OF THE COMMIT"`. 
  	* If you mistakely commit something, you can revert (undo) the last commit using `git revert HEAD`.
5. Once your project is complete, push the changes from your local clone to the origin (what you see online) using `git push origin NAME-OF-BRANCH`.
6. To merge your code, get online and submit a pull request from your branch in your fork to the main repo (see [this Git documentation for merging across forks](https://docs.github.com/en/github/collaborating-with-issues-and-pull-requests/creating-a-pull-request-from-a-fork))
7. To start a new project, go back to the `main` branch and start back at step 1

### Useful commands
Here is a simple list of useful commands when contributing with git:
* `git pull origin` gets updates from the online version of your fork
* `git checkout NAME-OF-BRANCH` changes to an existing branch
* `git checkout -b NAME-OF-BRANCH` creates a new branch and changes to it
* `git status` shows all files that have been changed or added to be committed
* `git add FILE OR FOLDER` stages changes to be committed
* `git commit -m "COMMENT TO DESCRIBE THE INTENTION OF THE COMMIT"` makes the change official on your local copy
* `git push origin NAME-OF-BRANCH` updates the the online version of your fork
* `git revert HEAD` undoes the last commit on a branch
* `git branch` lists the branches you have on your local copy; there will be a star by the branch you are currently on
* `git branch -d NAME-OF-BRANCH` will delete a given branch
* `git diff` shows the difference between local, unstaged changes and the official (committed) version
* `git checkout -- NAME-OF-FILE` deletes changes to a given file that have not been staged to commit
* `git checkout .` deletes all local changes in the repository that have not been added to the staging area
* `git clean -f` deletes untracked changes
* `git reset .` removes files from staging area before they have been committed
