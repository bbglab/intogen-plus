#!/bin/bash

# Check if the destination folder is provided as the last argument
if [ $# -lt 2 ]; then
  echo "Usage: $0 <source files>... <destination folder>"
  exit 1
fi

# Get the destination folder
destination="${@:$#}"

# Check if the destination folder exists
if [ ! -d "$destination" ]; then
  echo "Error: Destination folder does not exist"
  exit 1
fi

# Create symlinks for each source file
for source in "${@:1:$#-1}"; do
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
