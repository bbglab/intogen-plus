#!/bin/bash

# Function to create symlinks for files
create_symlinks_for_files() {
  local destination="$1"
  shift
  for source in "$@"; do
    # Check if the source file exists
    if [ ! -e "$source" ]; then
      echo "Error: $source does not exist"
      continue
    fi
    
    # Get the basename of the source file
    source_basename="$(basename "$source")"
    
    # Create the symlink
    ln -sfr "$source" "$destination/$source_basename"
    
    # Print a message indicating success
    echo "Created symlink for $source in $destination"
  done
}

# Function to create symlinks for folders
create_symlinks_for_folders() {
  local destination="$1"
  shift
  for source in "$@"; do
    # Check if the source folder exists
    if [ ! -d "$source" ]; then
      echo "Error: $source is not a valid folder"
      continue
    fi

    # Create symlinks for files inside the folder
    for file in "$source"/*; do
      create_symlinks_for_files "$destination" "$file"
    done
  done
}

# Main script
if [ $# -lt 2 ]; then
  echo "Usage: $0 <source files/folders>... <destination folder>"
  exit 1
fi

# Get the destination folder
destination="${@:$#}"

# Check if the destination folder exists
if [ ! -d "$destination" ]; then
  echo "Error: Destination folder does not exist"
  exit 1
fi

# Determine whether the input sources are files or folders
files=()
folders=()
for arg in "${@:1:$#-1}"; do
  if [ -f "$arg" ]; then
    files+=("$arg")
  elif [ -d "$arg" ]; then
    folders+=("$arg")
  else
    echo "Error: $arg is neither a file nor a folder"
  fi
done

# Create symlinks for files
create_symlinks_for_files "$destination" "${files[@]}"

# Create symlinks for files inside folders
create_symlinks_for_folders "$destination" "${folders[@]}"
