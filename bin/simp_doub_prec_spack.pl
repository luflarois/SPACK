#! /usr/bin/perl -w

#Ce script permet la transformation de la chaîne de charactère "E+" en
# "D+" dans *.f generes par SPACK (et E- en D-)

#----------------------------------------------------------------------------------------
#Création du tableau qui contient le nom des fichiers à modifier
#----------------------------------------------------------------------------------------
@fichs= ('./kinetic.f','./fexchem.f','fexprod.f','jacdchemdc.f','rates.f','fexloss.f');

$i=0;
#----------------------------------------------------------------------------------------

while ($i<@fichs) {

    $nomfic ="$fichs[$i]";
    $nomfic2="$fichs[$i]".".new";
    print "Ouverture du fichier "."$nomfic"."\n";
    open FIC , "$nomfic"    or die "Impossible d'ouvrir $nomfic\n";
    open FIC2, "+>$nomfic2" or die "Impossible d'ouvrir $nomfic\n";
    $i ++;
    while ($ligne=<FIC>) {
	print "$ligne";
	$ligne =~ s/E-/D-/g;
	$ligne =~ s/E\+/D\+/g;
	print FIC2 "$ligne";
    }
    rename "$nomfic2","$nomfic";
    close FIC;
    close FIC2;
}
