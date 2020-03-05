create table msa_2(orf_name varchar(255), species text, sequence longtext, _type text) charset utf8;
load data local infile 'SynORFan/msa_0214.csv' into table msa_2 fields terminated by ',' lines terminated by '\n';

alter table msa_2 add index(orf_name);


/*alter table synteny_long modify orf_name varchar(255);
alter table synteny_long add index(orf_name);*/

/*create table pairwise_best(orf_name varchar(255), species1 text, sequence1 longtext, species2 text, sequence2 longtext, type text, id int) charset utf8;
alter table pairwise_best add index(orf_name);
*/

