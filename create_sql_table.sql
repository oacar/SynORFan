## to save sql file use following on bash
mysqldump --single transaction -u... -p... db_name tb1 tb2 tb3 > flname.sql
create table msa_2(orf_name varchar(255), species text, sequence longtext, _type text) charset utf8;

load data local infile 'SynORFan/msa_0214.csv' into table msa_2 fields terminated by ',' lines terminated by '\n';

alter table msa_2 add index(orf_name);


/*alter table synteny_long modify orf_name varchar(255);
alter table synteny_long add index(orf_name);*/

/*create table pairwise_best(orf_name varchar(255), species1 text, sequence1 longtext, species2 text, sequence2 longtext, type text, id int) charset utf8;
alter table pairwise_best add index(orf_name);
*/
create table synteny_long(orf_name varchar(255), data varchar(255), Scer text, Spar text, Smik text, Skud text, Seub text, Suva text, Sarb text, Sjur text);
