#!/bin/bash

## how to use
### bash download_dx.sh $destination_local $target_dir_dx

export PATH="$HOME/.local/bin:$PATH"

destination_dir=$1
dx_project="project-GvFxJ08J95KXx97XFz8g2X2g"
dx_folder=$2

# Ensure destination directory exists
mkdir -p "$destination_dir"

# List files in the DNAnexus folder
for file in $(dx ls "$dx_project:$dx_folder"); do
    # Set the local file path
    local_file="$destination_dir/$file"

    # Check if the file already exists locally
    if [ -f "$local_file" ]; then
        echo "File exists locally: $local_file"

        # Get local file size
        local_size=$(stat -c%s "$local_file")

        # Get remote file size using dx stat
        remote_size=$(dx describe "$dx_project:$dx_folder/$file" --json | jq '.size')

        # Compare local and remote sizes
        if [ "$local_size" -eq "$remote_size" ]; then
            echo "Skipping $file (already complete)"
        else
            echo "File is incomplete. Resuming download..."
            dx download "$dx_project:$dx_folder/$file" -o "$local_file" -f
        fi
    else
        echo "Downloading $file..."
        dx download "$dx_project:$dx_folder/$file" -o "$local_file"
    fi
done


