mkdir -p plots.trash
mv plots/*.pdf plots.trash
set -e // Exit when any command fails
mkdir -p plots.an/ttbar/2016
rm -f plots.an/ttbar/2016/*.pdf
mkdir -p plots.an/ttbar/2017
rm -f plots.an/ttbar/2017/*.pdf
mkdir -p plots.an/ttbar/2018
rm -f plots.an/ttbar/2018/*.pdf
mkdir -p plots.an/ttbar/run2
rm -f plots.an/ttbar/run2/*.pdf
scons && ./run/higgsino/an_plot_syst_ttbar.exe --year 2016
mv plots/*.pdf plots.an/ttbar/2016
scons && ./run/higgsino/an_plot_syst_ttbar.exe --year 2017
mv plots/*.pdf plots.an/ttbar/2017
scons && ./run/higgsino/an_plot_syst_ttbar.exe --year 2018
mv plots/*.pdf plots.an/ttbar/2018
scons && ./run/higgsino/an_plot_syst_ttbar.exe --year 2016,2017,2018
mv plots/*.pdf plots.an/ttbar/run2

mkdir -p plots.an/ttbar.unblind/2016
rm -f plots.an/ttbar.unblind/2016/*.pdf
mkdir -p plots.an/ttbar.unblind/2017
rm -f plots.an/ttbar.unblind/2017/*.pdf
mkdir -p plots.an/ttbar.unblind/2018
rm -f plots.an/ttbar.unblind/2018/*.pdf
mkdir -p plots.an/ttbar.unblind/run2
rm -f plots.an/ttbar.unblind/run2/*.pdf
scons && ./run/higgsino/an_plot_syst_ttbar.exe --year 2016 --unblind
mv plots/*.pdf plots.an/ttbar.unblind/2016
scons && ./run/higgsino/an_plot_syst_ttbar.exe --year 2017 --unblind
mv plots/*.pdf plots.an/ttbar.unblind/2017
scons && ./run/higgsino/an_plot_syst_ttbar.exe --year 2018 --unblind
mv plots/*.pdf plots.an/ttbar.unblind/2018
scons && ./run/higgsino/an_plot_syst_ttbar.exe --year 2016,2017,2018 --unblind
mv plots/*.pdf plots.an/ttbar.unblind/run2
