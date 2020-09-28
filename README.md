# Drosophila genome analysis workflows
The scripts and Dockerfiles contained here are provided to enable the exact
reproduction of our computational pipelines.

## Downloading project files
The easiest way to access project files (for now) is to use Rclone, a command
line tool for managing cloud storage. Download it here:
https://rclone.org/downloads/

To set it up to work with our Box folder, copy the ```rclone``` binary into a 
directory in your ```$PATH```. Then, configure it by running ```rclone config```
and by choosing the following options. 
1. Add a new remote (```n```) 
1. Name it ```box```, or whatever you'd like
1. Choose option ```6``` for Box storage
1. Don't enter anything for Client Id
1. Don't enter anything for the Client Secret
1. Don't enter anything for the config.json location
1. Sub type is ```1``` or ```"user"```
1. Enter ```n``` and don't edit advanced config
1. If you are working on your own machine, enter ```y```, if you are trying to 
set up Rclone on a remote machine (e.g. the cluster), enter ```n``` 
   - If you're setting up on a local machine, an authorization web page should
     automatically appear. If it does not, go to the address provided by Rclone.
   - If you're setting up on a remote machine, you will need to have a copy of
     ```rclone``` installed on your local machine. At this step you will need to
     run the command ```rclone authorize "box"``` on your local machine, 
     authorize using the web page that shows up, then provide the authorization
     key to the remote machine.

Once ```rclone``` has been set up in this manner, downloading a Box directory is
as easy as running:

```bash
rclone copy -P box:100x100/assemblies/repeat_masked/fastas /path/to/local/dir
```

If you have download some files already and want to skip them:
```bash
rclone sync -P box:100x100/assemblies/repeat_masked/fastas /path/to/local/dir
```

Full list of commands: https://rclone.org/commands/