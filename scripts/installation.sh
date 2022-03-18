#!/bin/sh

cd $HOME

BASH_RC="$HOME"/.bashrc
BASH_FILE="$HOME"/.profile
BASH_PROFILE="$HOME"/.bash_profile

if [ -f "$BASH_RC" ]; then
    printf "\\n"
    printf "Initializing PyOrb in %s\\n" "$BASH_RC"
    printf "\\n"

cat <<EOF >> "$BASH_RC"
# added by PyOrb installer
export PYORBPATH="$HOME/pyorb"
export PATH="\$PYORBPATH/bin:\$PATH"
# added by PyOrb installer
EOF

source "$HOME"/.bashrc

elif [ -f "$BASH_FILE" ]; then
    printf "\\n"
    printf "Initializing PyOrb in %s\\n" "$BASH_FILE"
    printf "\\n"

cat <<EOF >> "BASH_FILE"
# added by PyOrb installer
export PYORBPATH="$HOME/pyorb"
export PATH="\$PYORBPATH/bin:\$PATH"
# added by PyOrb installer
EOF

source "$HOME"/.profile

else
    printf "\\n"
    printf "Initializing PyOrb in %s\\n" "$BASH_PROFILE"
    printf "\\n"

cat <<EOF >> "$BASH_PROFILE"
# added by PyOrb installer
export PYORBPATH="$HOME/pyorb"
export PATH="\$PYORBPATH/bin:\$PATH"
# added by PyOrb installer
EOF

source "$HOME"/.bash_profile
fi


pip install numpy
pip insatll matplotlib
pip install pandas
pip install seaborn


curl -L -o aichem.zip  https://github.com/sunxb05/aichem/zipball/master/
PYORBPREZIP="$HOME/aichem.zip"
unzip "$PYORBPREZIP"
rm -rf "$PYORBPREZIP"
mkdir $HOME/pyorb
mv $HOME/sunxb05-aichem*/pyorb $HOME/pyorb
mv $HOME/sunxb05-aichem*/bin $HOME/pyorb
rm -r $HOME/sunxb05-aichem*
