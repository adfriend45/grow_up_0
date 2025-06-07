cd $HOME/grow_up_0/src
make
cd $HOME/grow_up_0/scripts
mv $HOME/grow_up_0/src/tgu.exe .
nice -19 ./tgu.exe
cd $HOME/grow_up_0/data
python3 plot_output.py &
