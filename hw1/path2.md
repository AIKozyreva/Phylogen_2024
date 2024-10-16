#### Добаавляем материалы для Неандертальцев и Денисовцев, переделываем выравнивание и дерево. 

#### Потом создаём метаданные нужного формата и подгружаем всю красоту в ITOL.
```
План был хорош, но есть нюанс - у меня никогда ничего не идёт по плану. 
Начнём с создания выравниваний для материалов Неандертальцев и Денисовцев, а затем объеденим это с выравниванием, по которому строилось дерево в path1. 
```
__________________________________
## Наскребаем метаданные для дерева из NCBI gbk файлов. 
В прошлой части мы слили вместе разные фасты из какого-то датасета (мы не знаем откуда они, но надеемся, что они все имеют стандартный фаста-заголовок с ncbi, у такого после знака `>` стоит ncbi-accession уникальный набор букв и числел). Мы эти все фасты сливали вместе, префращали в мультифасту, выравнивали записи в ней, получали таким образом мультивыравнивание, а потм по нему строили дерево и имели на ветках уникальные наборы цифро-букв, а именно ncbi-accessions. Чтобы поместить на дерево не цифро-буквы, а какие-то адекватные слова, имеющие значения, и тем самым сделать дерево читабельным - нам нужны метаданные о фаста-записях, которые мы использовали для построения выравнивания. Метаданные хранятся в базе ncbi и могут быть представлены в куче форматов. Иногда загружающие последовательность люди заполняют поля отдельно в установленном формате, тогда вы можете так или иначе попросить у базы ncbi только эти данные, а иногда данные нужно искать, потому что это нечто особенное (оно вносится в поле "note" или другие, например раздела SOURCE, в свободном формате). Для особенной инфы нам нужно просмотреть глазами часто GBK файл (это то, что вы видите открывая какую-то запись на сайте ncbi). 

![image](https://github.com/user-attachments/assets/4a947641-9a55-422a-8606-1074288af2b1)


Для большого кол-ва образцов никто, конечно, смотреть GBK не будет - их массово скачивают, а потом парсят, пользуясь какими-либо текстовыми паттернами и/или регулярными выражениями, пытаясь достать нужную инфу. Вот именно это мы и сделаем в этом разделе через боль слёзы и я так ненавижу питон прости гсосподи, ну ладно R хуже. 

**Шаг 1 Как достать accessions NCBI из мультифасты**
Эти аксешн нуджны, зачем будет ясно в шаге2, но где то надо иметь их списком, поэтому будем их доставать. Открываем Python, загружаем в память вашу мультифасту или просто протягиваемся до неё (у мменя ноутбук питона на серве поднят, а файлы лежали там же), и пытаемся извлечь максимум пользы из заголовков каждой фаста-записи в мультифасте. 

```
import pandas as pd
from Bio import SeqIO

fasta_file = '/mnt/230924_phylo_hw/comb_part2_NDwmtHT.fasta'
data = []
for record in SeqIO.parse(fasta_file, "fasta"):
    header = record.description
    accession = header.split('Homo')[0]    
    if "haplogroup" in header:
        haplogroup = header.split("haplogroup")[1].split('mitochondrion')[0].strip()
    else:
        try:
            if "sapiens" in header:
                haplogroup = header.split('sapiens')[1].split('mitochondrion')[0].strip()
            else: 
                haplogroup = header.split('Homo')[1].split('from')[1].strip()
        except IndexError:
            haplogroup = "unknown"            
    data.append([accession, haplogroup])

haplo_meta_df = pd.DataFrame(data, columns=["accession", "haplogroup"])
haplo_meta_df.to_csv('/mnt/230924_phylo_hw/061024_part2/haplo_meta_df.tsv',sep="\t")
```
Всё - сохранили метаданные (ну как уж получилось) в табличку. Посмотрим чего там насохранялось. 
![image](https://github.com/user-attachments/assets/287f5ab1-6e26-47bc-bd9f-355ec2d234f0) 

Угу. Это безобразие результат того, что заголовки фасты выглядели по-разному, но имели какие-то общие черты (у всех строго в начале был ncbi-accession, а 
 дальше в свободном формате). Это плохо, потому что становится неудобно парсить, и вы получаете мою табличку. Не страшно в нашем случае, тк это подтверждает необходимость доставать метаданные не из фаста файла, а сложным путём через запрос сведений в бд ncbi для каждой уникальной фасты. 
 
   важнейший принцип биоинформатики be like: **Можно безобразно - главное однообразно** :) 

Мы получили хорошую колонку Accessions - по факту это как раз то, что нужно -список - чтоб по этим id завпросить у бд gbk файлы и их прочитать. (в моём случае я их ещё сохраню для спокойствия души, но этого делать НЕ стоит, если у вас мало места). 


**Шаг 2 Как достать метаданные с NCBI**
Теория: у NCBI есть два варика программного доступа для массового или автоматического скачивания данных. 

Вариант 1: Это либо API (программный интерфейс) -> оно называется E-utils, оно работает как часть библиотеки Biopython, и я уже даже понимаю как оно работает (примерно), но лучше от этого не становится.

Вариант 2: Это пакет программ "ncbi-datasets-cli". Ставится через конду `conda install -c conda-forge ncbi-datasets-cli` и там есть программки datasets и dataformat - замечательные вещи, через которые удобно запрашивать большое кол-во типов данных и большое по объёму количество, НОО но только из баз genomes, gene, viruses NCBI. То есть доступ к случайной записи, которую вы через веб интерфейс находите в базе NCBI nucleotides (и там лежит метаинфа и сама например запись с частью какого-т гена, то как бы нет, он просто не пропустит запрос - выдаст ошибку bad request с кодом 400 - потому что такой базы он не знает).

мб недостаток Варианта2 - это специально, а мб они когда-то допилят доступ с этих тулов ко всей своей библиотеке, но пока для доступа к любому отделу NCBI, который не genome/gene/viruses - надо работать через E-utils.

Документация на E-utils NCBI: https://www.ncbi.nlm.nih.gov/books/NBK25498/#chapter3.ESearch__ESummaryEFetch (смотреть примеры использования Efetch Esummary); https://biopython.org/docs/1.76/api/Bio.Entrez.html 

Документация на ncbi-datasets-cli: https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/ 

Факт: для того, чтобы достать любые данные через любой из вариантов выше вам нужно указать какой-то идентификатор записи, что прога вообще должна спросить у базы. Чем больше параметров вы указываете, тем уже будет вопрос проги к базе - всё как если б вы выставляли фильтры ручками на сайте ncbi, но только происходит под капотом. В качестве таких идентификаторов могут выступать или accessions или ncbi taxon id. Обычно именно эти значения. У нас как раз accessions теперь есть из предыдущего шага (колонка в датафрейме).

В коде ниже облажалась с парсингом гбк файла, надо будет переписать когда-нибудь, потому что я напрасно ориентировалась только на поле note, надо было ещё geo_loc добавить
```
from Bio import Entrez, SeqIO
import os
import pandas as pd
import re
import os

Entrez.email = "anfisakozyr@gmail.com"
Entrez.api_key = "262c18e128672487a2ad30be17815ca37008"

input_tsv = "/mnt/230924_phylo_hw/061024_part2/haplo_meta_df.tsv"
output_dir = "/mnt/230924_phylo_hw/061024_part2/gbk_files"
output_metadata_file = "/mnt/230924_phylo_hw/061024_part2/metadata.tsv"

# Функция для Esearch 
def search_accession(accession):
    search_handle = Entrez.esearch(db="nucleotide", term=accession, usehistory="y")
    search_results = Entrez.read(search_handle)
    search_handle.close()
    print (f" \n search result for {accession} is \n {search_results}")
    return search_results

# Функция для EFetch
def fetch_gbk_file(accession, webenv, query_key):
    efetch_handle = Entrez.efetch(db="nucleotide", query_key=query_key, WebEnv=webenv, rettype="gb", retmode="text")
    gbk_filename = os.path.join(output_dir, f"{accession}.gbk")
    print(efetch_handle.read())
    with open(gbk_filename, "w") as gbk_file:
        gbk_file.write(efetch_handle.read())
    efetch_handle.close()
    return gbk_filename

def parse_gbk_file(gbk_filename):
    handle = open(gbk_filename, "r")
    record = SeqIO.read(handle, "genbank")
    handle.close()
    return record

def extract_note_info(note):
    country = ethnicity = "unknown"
    for eth in ethnicities:
        if eth.lower() in note.lower():
            ethnicity = eth
            break
    for country_name in countries:
        if country_name.lower() in note.lower():
            country = country_name
            break
    return ethnicity, country

os.makedirs(output_dir, exist_ok=True)
df = pd.read_csv(input_tsv, sep='\t', header=None, names=["accession", "haplogroup"])
accessions = df['accession']
metadata_rows = []

for accession in accessions:
    try:
        print(f"Searching for accession: {accession}")
        
        # Выполняем поиск - esearch
        search_results = search_accession(accession)
        if search_results["Count"] == "0":
            print(f"No records found for {accession}")
            continue
        webenv = search_results["WebEnv"]
        query_key = search_results["QueryKey"]

        # Выпаолняем запрос efetch - получаем доступ и сохраняем GBK
        print(f"Fetching GBK file for {accession}")
        gbk_filename = fetch_gbk_file(accession, webenv, query_key)
        
        # Парсинг GBK файла
        record = parse_gbk_file(gbk_filename)
        for feature in record.features:
            if feature.type == "source":
                haplogroup = feature.qualifiers.get("haplogroup", ["unknown"])[0]
                organelle = feature.qualifiers.get("organelle", ["unknown"])[0]
                note = feature.qualifiers.get("note", ["unknown"])[0]

                # Извлекаем страну и национальность из note
                ethnicity, country = extract_note_info(note)
                note_info = f"Ethnicity: {ethnicity}, Country: {country}"
                break

        # Добавление метаданных в список
        metadata_rows.append([accession, haplogroup, organelle, note_info])
    except Exception as e:
        print(f"Error fetching data for {accession}: {e}")
        metadata_rows.append([accession, "error", "error", str(e)])

# Сохранение метаданных
metadata_df = pd.DataFrame(metadata_rows, columns=["accession", "haplogroup", "organelle", "note_info"])
metadata_df.to_csv(output_metadata_file, sep='\t', index=False)

print(f"Metadata saved to {output_metadata_file}")
print(f"GBK files saved to {output_dir}")

```
