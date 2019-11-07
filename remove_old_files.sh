find -name "*subalignment*" | while read line; do rm $line; done;
find -name "*data*" | while read line; do rm $line; done;
find -name "*AATranslation*" | while read line; do rm $line; done;
find -name "*muscle*" | while read line; do rm $line; done;

find -name "S*" | while read line; do rm -r $line; done;
