## prerekvizity

technické: python, pip, virtualenv
iné: zlib1g-dev

## Príprava

Je potrebné vytvoriť virtualenv a rozbehnúť potrebné knižnice.

0.  $ virtualenv venv

1.  $ source venv/bin/activate

2.  (venv) $ pip install -r requirements.txt

## Použitie

    (venv) $ python3 refmatch.py -r test/reference.fa -s test/

Tento príkaz spustí program refmatch a referenciu zo súboru ./test/reference.fa vyhľadá vo .fast5 súboroch
z priečinka test(takisto aj v jeho podpriečinkoch)

## Výstup

Program vypíše výstup do súboru ./out.txt . Výstup pozostáva z niekoľkých riadkov, na každom sa nachádza
cesta k .fast5 súboru, v ktorom bola identifikovaná zhoda. Potom na tomto riadku nasleduje niekoľko
čísel, ktoré symbolizujú signál, ktorý zodpovedá našej referencii v tomto .fast5 súbore.
