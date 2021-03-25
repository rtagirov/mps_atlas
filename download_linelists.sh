#!/bin/bash

mkdir linelists
cd linelists || exit
wget -q --content-disposition "https://owncloud.gwdg.de/index.php/s/lnC5Sjs9K0U7VEt/download?path=%2F&files=diatomicsdf.bin&downloadStartSecret=il4b99hwr3n" -O diatomicsdf.bin
wget -q --content-disposition "https://owncloud.gwdg.de/index.php/s/lnC5Sjs9K0U7VEt/download?path=%2F&files=highlinesdf.bin&downloadStartSecret=76glg6xrjmo" -O highlinesdf.bin
wget -q --content-disposition "https://owncloud.gwdg.de/index.php/s/lnC5Sjs9K0U7VEt/download?path=%2F&files=tiolinesdf.bin&downloadStartSecret=4xtlwyc46gi" -O tiolinesdf.bin
wget -q --content-disposition "https://owncloud.gwdg.de/index.php/s/lnC5Sjs9K0U7VEt/download?path=%2F&files=nltelinesdf.bin&downloadStartSecret=l2pfml7rsd" -O nltelinesdf.bin
wget -q --content-disposition "https://owncloud.gwdg.de/index.php/s/H69BNWi7jI2jQGY/download" -O h2olinesdf.bin
wget -q --content-disposition "https://owncloud.gwdg.de/index.php/s/ZrytLF5MZKqiFsk/download" -O lowlinesdf.bin
