**Статья на воспроизведение: Francis WR, Shaner NC, Christianson LM, Powers ML, Haddock SH. Occurrence of Isopenicillin-N-Synthase Homologs in Bioluminescent Ctenophores and Implications for Coelenterazine Biosynthesis. PLoS One. 2015;10(6):e0128742. Published 2015 Jun 30. doi:10.1371/journal.pone.0128742**

**Данные для анализа:** Последовательности из статьи, включающие для разных таксонов ряд транскриптов, потенциально кодирующих интересующие авторов пептиды с FYY на конце и\или потенциально участвующие в цикле синтеза коцелентироза, который является не окисленным вариантом субстрата для ферментативной реакции люцефирин-люциферазного взаимодействия, которая в разных морских организмах приводит к биолюминесценции. "_Ctenophore sequences used in analysis can be found at GenBank, with accessions: KM233765-KM233833_".

""_Здесь мы искали гены, кодирующие «FYY» из транскриптомов светящихся гребневиков. Нас также интересовали гены, которые потенциально могли бы выполнять шаги циклизации, обсуждавшиеся выше. Мы идентифицировали гены-кандидаты, которые присутствовали в транскриптомах светящихся видов и отсутствовали для несветящихся видов. Мы сравниваем эти белки с белками из геномов родственных животных и показываем, что эта группа белков высококонсервативна даже среди отдаленно родственных гребневиков, что ожидается для критических биологических процессов._""

![image](https://github.com/user-attachments/assets/a0aa7377-726c-4a31-92c4-a5029eb78ebb)

изначально группа анализировала полные транскриптомы, но так как это входит в смысл моей дипломной, которой ещё пока нет, а не этого проектка, то мне нет смысла повторять самуую сложную часть. Повторим только филогенетику. 

**Аутгруппы:** Сложно, так как они строили дерево фактически по последовательностям различных белков для разных гребневиков. Я включила также ген, кодирующий изопенницилин N синтетазу в бактериях, но вряд ли это уместно, так как они на одном из этапов искали "на что похожи обнаруженные fyy последовательности", и совпадение BLAST было как раз на бактериальную изопенницилин N-синтетазу. то есть, что в общем не новость, белковые элементы наиболее консервативны для единых функций, даже в оч удалённых таксонах. А брать что угодно, кодирующее совсем левый белок - ну как-то странно, очевидно, что будет просто не похож, он же даже не выровняется. В итоге добавила в третье дерево гомолог FYY-пептидов из бактерии (что далеко от гребневиков): 1) >NZ_AHAM01000058.1:Mesorhizobium alhagi CCNWXJ12-2, isopenicillin N synthase. 

Вот так выглядит выравнивание нуклеотидных последовательностей из N гребневиков, у которых были обнаружены FYY-концевые участки. Выравнивание построено без концевых групп и без удаления gap, обрезано только по краям. Не очень хорошее выравнвивание объяснимо - пусть это и мотивы FYY на конце, но в этой же статье пишут, что они относятся к разным классам белоков, участвующих в функциях похожих, но разные классы белков - это сильно большая разница. 

![image](https://github.com/user-attachments/assets/6a8530d2-2ee5-41a1-ac5a-aa5afb5ea78e)

Дерево в статье, демонстрирующее сродство FYY-последовательностей вне зависимости от таксона. То есть они исходно искали все воследовательности с мотивом FYY, потом посмотрели, что часть таких являются гомологами изопенницилин-н-синтетаз и ещё одного класса белков. Потом взяли все последовательности, похожие на эти три группы в транскриптоме каждого исследуемого таксона. Дерево внизу (статейное) - демонстрирует результаты кластеризации по подобию - вне зависимости от таксона последовательност разделяются по типам белков, гомологами которых они являются.

![image](https://github.com/user-attachments/assets/1aff13f6-96b4-4b03-a9ff-9a770c89cfe6)


![image](https://github.com/user-attachments/assets/0b4116a2-d98f-4981-9619-02dd580f20e1)


Я сделала 3 дерева: 1) полностью на сыром выравнивании всех участвовавших в статье последовательностей-FYY, методом правдоподобия с моделью TVM+F+R5. 2) Для группы гомологов isopnicillin-n-synthase бактерий из разных изучаемых таксонов дерево. 3) Для группы гомологов эукариотической 2OGFe1(2)-oxygenase из разных изучаемых таксонов дерево.

Команды аналогичные для каждого этапа. В общем виде показаны ниже. 
```
mafft --kimura 1 --clustalout --maxiterate 1000 --thread 4 --auto file.fasta > alignment_file.aln
iqtree2 -s alignment_file.aln -bb 1000 -m MFP -nt 3

```
______________________________________________________________________________________________________________________________________________________________

#### 1) На выравнивании всех участвовавших в статье полследовательностей FYY, методом правдоподобия с моделью TVM+F+R5.

![OBJBLUFrQBRUrRLMZq1WkA](https://github.com/user-attachments/assets/0e87d8e1-82c0-42ed-8685-674019bd902e)

MLXXX = Mnemiopsis leidyi - model organism's FYY-peptides, which have have had different classes hits among their genome data. Hormiphora californensis - не светится.

______________________________________________________________________________________________________________________________________________________________

#### 2) Дерево для последовательностей-FYY, гомологичных isopenicillin-N-synthase бактерий, в разных гребневиках.

Дерево представлено на рисунке ниже. Дерево с подкреплённым к нему выравниванием представлено в файле "FYY-pepides_alignment_phylo.pdf"
Дерево с подбором эволюционной модели построено методом правдоподобия с моделью TPM2+F+R4. 

![Eg5NvmqY-s-iuQ-9-6nFRQ](https://github.com/user-attachments/assets/3de561e3-6a1d-4c8f-9866-bcf7877656aa)

______________________________________________________________________________________________________________________________________________________________

#### 3) Дерево для последовательностей-FYY, гомологичных 2OGFe1(2), в разных гребневиках.

Дерево представлено на рисунке ниже. Дерево с подкреплённым к нему выравниванием представлено в файле "2OGFe12-pepides_alignment_phylo.pdf"
Дерево с подбором эволюционной модели построено методом правдоподобия с моделью TVMe+R4. 


![image](https://github.com/user-attachments/assets/903c7b99-0edd-4cc5-b1be-335b8c493ccd)

______________________________________________________________________________________________________________________________________________________________

![image](https://github.com/user-attachments/assets/93c82ef0-82f3-461c-b3c9-ce3c143ff069)


"_Interestingly, nearly all of the top hits for all of the proteins were to a 2OG-Fe(II) oxygenase from the ciliate Oxytricha trifallax. This was surprising since ciliates are unicellular eukaryotes and are not closely related to ctenophores. In a more restricted search using the Uniprot/Swissprot database, the top BLAST hits for many of the FYY proteins were to the same set of isopenicillin-N-synthase (IPNS) homologs, mostly from bacteria. These proteins are members of a group of Fe-dependent oxygenases that include IPNS and deacetoxycephalosporin C synthase (DAOCS). These are the enzymes responsible for the heterocycle-forming steps of penicillin biosynthesis and the ring expansion in cephalosporin biosynthesis, respectively, and therefore were considered even stronger candidates for involvement in cyclization of FYY to coelenterazine_".

```
Species	Luminous? Y/N	
Bathocyroe fosteri	Yes	
Bathyctena chuni	Yes	
Beroe abyssicola	Yes	
Beroe forskalii	Yes
Bolinopsis infundibulum	Yes	
Charistephane fugiens	Yes
Dryodora glandiformis	Yes
Euplokamis dunlapae	Yes
Haeckelia rubra	Yes	
Hormiphora californensis	No
Lampea lactea	Yes	
Lampocteis cruentiventer	Yes	
Ocyropsis maculata	Yes	
Thalassocalyce inconstans	Yes
Undescribed ctenophore B	Yes
Undescribed ctenophore C	Yes	
Undescribed ctenophore N1	Yes	
Undescribed ctenophore N2	Yes	
Undescribed ctenophore T	Yes	
Undescribed ctenophore V	Yes	
Undescribed ctenophore W	Yes	
Velamen parallelum	Yes	
```



P.S. Про программы. 

Вравнивание mafft (визуализация - ugene, хотя учитывая небольшой размер последовательностей можно было всё сразу в ugene делать, тк никаких особенных параметров я не ставила). Дерево строила с использованием iqtree2, потому что я уже умела его использовать (+-), а ещё у него еть ModelFinder, это модуль, который подбирает подходящую модель эволюции и оптимизирует под неё параметры при выстраивании дерева, чтобы не надо было доп.тупы запускать. Визуализация дерева ITOL. Подготовка файлом с метаданными для декорирования дерева - python.
