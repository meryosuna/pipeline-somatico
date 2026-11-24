# Pipeline somático com GATK4 Mutect2

Atividade prática – Bioinformática Aplicada à Genômica Médica
Hospital Israelita Albert Einstein

## 1. Objetivo

Compreender as principais etapas de um pipeline de análise de variantes somáticas utilizando o GATK4 Mutect2, desde a preparação e indexação da referência até a chamada e filtragem das variantes.

## 2. Ferramentas utilizadas

- GATK 4.2.2.0
- SAMtools
- Bash
- GitHub Codespaces
- Referência hg19
- gnomAD
- Panel of Normals (PoN)

## 3. Preparação da referência

### 3.1. Download dos cromossomos 13 e 19

```bash
wget -c https://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/chr19.fa.gz
wget -c https://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/chr13.fa.gz
```

Conferir se os arquivos foram salvos:

```bash
ls -lh chr13.fa.gz chr19.fa.gz
```

### 3.2. Montagem da referência hg19

```bash
zcat chr13.fa.gz chr19.fa.gz | sed 's/^>chr/>/' > hg19.fa
```

Conferir os cabeçalhos:

```bash
grep '^>' hg19.fa
```

### 3.3. Indexação da referência

O SAMtools permite criar um índice para facilitar o acesso às sequências da referência.

```bash
samtools faidx hg19.fa
cat hg19.fa.fai
```

### 3.4. Criação do dicionário da referência

```bash
./gatk-4.2.2.0/gatk CreateSequenceDictionary \
  -R hg19.fa \
  -O hg19.dict
```

Conferir as sequências:

```bash
grep '^@SQ' hg19.dict
```

### 3.5. Geração da lista de intervalos

No procedimento original, os intervalos são definidos utilizando o ScatterIntervalsByNs.

```bash
./gatk-4.2.2.0/gatk ScatterIntervalsByNs \
  -R hg19.fa \
  -O hg19.interval_list \
  -OT ACGT
```

```bash
awk '!/^@/ {print $1}' hg19.interval_list | sort -n | uniq -c
```

Durante a execução da atividade, foi utilizado um arquivo BED como alternativa à lista de intervalos.

```bash
printf "13\t0\t115169878\n19\t0\t59128983\n" > hg19_regions.bed
```

## 4. Análise tumor-normal: WP190/WP191

### 4.1. Conferência dos arquivos e identificação das amostras

```bash
ls -lh *wp190* *wp191* af-only-gnomad-chr13-chr19.vcf.gz*

samtools view -H tumor_wp190.bam | grep '^@RG'
samtools view -H tumor_wp191.bam | grep '^@RG'
```

### 4.2. Chamada de variantes com Mutect2

O Mutect2 permite identificar variantes somáticas utilizando a amostra tumoral e a amostra normal pareada.

```bash
./gatk-4.2.2.0/gatk Mutect2 \
  -R hg19.fa \
  -I tumor_wp190.bam \
  -I tumor_wp191.bam \
  -normal WP191 \
  --germline-resource af-only-gnomad-chr13-chr19.vcf.gz \
  -L hg19_regions.bed \
  -O WP190.somatic.vcf.gz
```

### 4.3. Geração das tabelas de contagem

```bash
./gatk-4.2.2.0/gatk GetPileupSummaries \
  -I tumor_wp190.bam \
  -V af-only-gnomad-chr13-chr19.vcf.gz \
  -L hg19_regions.bed \
  -O WP190.table
```

```bash
./gatk-4.2.2.0/gatk GetPileupSummaries \
  -I tumor_wp191.bam \
  -V af-only-gnomad-chr13-chr19.vcf.gz \
  -L hg19_regions.bed \
  -O WP191.table
```

```bash
head -n 5 WP190.table WP191.table
```

### 4.4. Estimativa de contaminação

```bash
./gatk-4.2.2.0/gatk CalculateContamination \
  -I WP190.table \
  -matched WP191.table \
  -O WP190.contamination.table
```

```bash
cat WP190.contamination.table
```

### 4.5. Filtragem das variantes

```bash
./gatk-4.2.2.0/gatk FilterMutectCalls \
  -R hg19.fa \
  -V WP190.somatic.vcf.gz \
  --contamination-table WP190.contamination.table \
  -O WP190.filtered.vcf.gz
```

Contagem das variantes de acordo com a coluna FILTER:

```bash
zcat WP190.filtered.vcf.gz | awk '!/^#/ {print $7}' | sort | uniq -c
```

## 5. Análise tumor-normal: WP017/WP018

### 5.1. Conferência das amostras

```bash
ls -lh *wp017* *wp018*

samtools view -H tumor_wp017.bam | grep '^@RG'
samtools view -H tumor_wp018.bam | grep '^@RG'
```

### 5.2. Chamada de variantes

```bash
./gatk-4.2.2.0/gatk Mutect2 \
  -R hg19.fa \
  -I tumor_wp017.bam \
  -I tumor_wp018.bam \
  -normal WP018 \
  --germline-resource af-only-gnomad-chr13-chr19.vcf.gz \
  -L hg19_regions.bed \
  -O WP017.somatic.vcf.gz
```

### 5.3. Geração das tabelas de contagem

```bash
./gatk-4.2.2.0/gatk GetPileupSummaries \
  -I tumor_wp017.bam \
  -V af-only-gnomad-chr13-chr19.vcf.gz \
  -L hg19_regions.bed \
  -O WP017.table
```

```bash
./gatk-4.2.2.0/gatk GetPileupSummaries \
  -I tumor_wp018.bam \
  -V af-only-gnomad-chr13-chr19.vcf.gz \
  -L hg19_regions.bed \
  -O WP018.table
```

```bash
head -n 5 WP017.table WP018.table
```

### 5.4. Estimativa de contaminação

```bash
./gatk-4.2.2.0/gatk CalculateContamination \
  -I WP017.table \
  -matched WP018.table \
  -O WP017.contamination.table
```

```bash
cat WP017.contamination.table
```

### 5.5. Filtragem das variantes

```bash
./gatk-4.2.2.0/gatk FilterMutectCalls \
  -R hg19.fa \
  -V WP017.somatic.vcf.gz \
  --contamination-table WP017.contamination.table \
  -O WP017.filtered.vcf.gz
```

```bash
zcat WP017.filtered.vcf.gz | awk '!/^#/ {print $7}' | sort | uniq -c
```

## 6. Identificação das variantes PASS

As variantes classificadas como PASS são aquelas que passaram pelos filtros aplicados pelo GATK.

```bash
for vcf in WP190.filtered.vcf.gz WP017.filtered.vcf.gz; do
  printf '\n%s\n' "$vcf"
  zcat "$vcf" | awk '/^#CHROM/ || (!/^#/ && $7 == "PASS")'
done
```

## 7. Análise utilizando Panel of Normals (PoN)

O Panel of Normals auxilia na identificação de artefatos recorrentes. Nesta estratégia, não foi utilizada a amostra normal pareada.

Exemplo da chamada de variantes para WP190:

```bash
./gatk-4.2.2.0/gatk Mutect2 \
  -R hg19.fa \
  -I tumor_wp190.bam \
  --germline-resource af-only-gnomad-chr13-chr19.vcf.gz \
  --panel-of-normals Mutect2-exome-panel.vcf \
  -L hg19_regions.bed \
  -O WP190.pon.somatic.vcf.gz
```

### 7.1. Estimativa de contaminação

```bash
./gatk-4.2.2.0/gatk GetPileupSummaries \
  -I tumor_wp190.bam \
  -V af-only-gnomad-chr13-chr19.vcf.gz \
  -L hg19_regions.bed \
  -O WP190.pon.table
```

```bash
./gatk-4.2.2.0/gatk CalculateContamination \
  -I WP190.pon.table \
  -O WP190.pon.contamination.table
```

### 7.2. Filtragem

```bash
./gatk-4.2.2.0/gatk FilterMutectCalls \
  -R hg19.fa \
  -V WP190.pon.somatic.vcf.gz \
  --contamination-table WP190.pon.contamination.table \
  -O WP190.pon.filtered.vcf.gz
```

O mesmo procedimento foi aplicado à amostra WP017, utilizando os respectivos arquivos e identificadores.

## 8. Análise da amostra JAK2 (WP043/WP044)

Como o gene JAK2 está localizado no cromossomo 9, foi necessário preparar uma referência específica para esse cromossomo.

### 8.1. Preparação e indexação

```bash
wget -c https://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/chr9.fa.gz

gunzip -c chr9.fa.gz | sed '1s/^>chr9/>9/' > hg19_chr9.fa

samtools faidx hg19_chr9.fa

./gatk-4.2.2.0/gatk CreateSequenceDictionary \
  -R hg19_chr9.fa
```

```bash
printf "9\t0\t141213431\n" > hg19_chr9_regions.bed
```

### 8.2. Adaptação do script

O script original foi adaptado para utilizar a referência, os intervalos e o recurso gnomAD do cromossomo 9.

```bash
sed \
-e 's/af-only-gnomad-chr13-chr19.vcf.gz/af-only-gnomad-chr9.vcf.gz/' \
-e 's/hg19_regions.bed/hg19_chr9_regions.bed/' \
-e 's/hg19.fa/hg19_chr9.fa/' \
run.chr13-chr19.sh > run.chr9.sh
```

Execução do pipeline tumor-normal:

```bash
sh run.chr9.sh tumor_JAK2.bam normal_JAK2.bam WP043 WP044
```

Contagem das variantes PASS:

```bash
zcat WP043.filtered.vcf.gz | awk '!/^#/ && $7=="PASS"' | wc -l
```

## 9. Resultados

### Comparação entre as estratégias

| Amostra tumoral | Tumor-normal (PASS) | PoN (PASS) |
|---|---:|---:|
| WP190 | 2 | 3 |
| WP017 | 1 | 1 |
| WP043 (JAK2) | 1 | Não concluído |

### Variantes identificadas nas amostras WP190 e WP017

| Amostra | Cromossomo | Posição | Tipo |
|---|---:|---:|---|
| WP190 | 13 | 28608257 | Inserção |
| WP190 | 13 | 28624322 | SNV |
| WP017 | 19 | 13054564 | Deleção |

Na análise com PoN, foi identificada uma variante adicional na amostra WP190, na posição 28602226 do cromossomo 13.

Para verificar essa variante, foi utilizado o seguinte comando:

```bash
samtools mpileup \
  -f hg19.fa \
  -r 13:28602226-28602228 \
  tumor_wp190.bam tumor_wp191.bam
```

Foi observada evidência da deleção tanto na amostra tumoral quanto na normal, sugerindo uma possível origem germinativa.

## 10. Conclusão

A atividade permitiu compreender as principais etapas do pipeline de análise de variantes somáticas, incluindo a preparação da referência, a chamada de variantes, a estimativa de contaminação e a filtragem do VCF.

A comparação entre as estratégias tumor-normal e PoN demonstrou a importância da amostra normal pareada na diferenciação entre variantes somáticas e germinativas.

Além disso, a classificação PASS indica que a variante passou pelos filtros técnicos, sendo necessária uma avaliação adicional para determinar sua relevância clínica.

## Referências

- Material disponibilizado em aula: https://github.com/renatopuga/somatico
- GATK: https://gatk.broadinstitute.org/
- UCSC Genome Browser: https://genome.ucsc.edu/
- gnomAD: https://gnomad.broadinstitute.org/
