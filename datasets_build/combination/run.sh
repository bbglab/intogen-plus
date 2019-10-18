# Define your key like this:
export COSMIC_KEY=$(echo "email@example.com:mycosmicpassword" | base64)

./run_negative_set.sh
./run_cgc.sh