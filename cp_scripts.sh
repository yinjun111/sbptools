#synchronize scripts in the dev codes folder with the running folder
rsync -a -d --delete /home/jyin/Projects/Pipeline/sbptools/ /apps/sbptools/

#make everything execuable +x
chmod 755 -R /apps/sbptools/
