Pokretanje skripti:
1. Capturing
2. Unpacking
3. Processing

Ono na šta trebaš paziti je raspored foldera.
Nisam nigdje napravio to da automatski radi foldere, pa se to mora ruèno podesiti. 
U tree.txt imas prikazano cijelo stablo. Znam da bi bilo lakše da je to automatski,
ali sad nemam vremena za to. Još samo napomena, to stablo ide od .\data\Scenes\

Kao što je najavljeno, MPS slike se spremaju u folder, ostale idu u èelije ili matrice.
Rezultati se spremaju isto u èeliju. 
Nisam probao skriptu 'Processing' za sve slike, samo za jednu. Sutra æu to izvrtiti.

Devignetting još nije popravljen. 
Funkcija 1: Cut_Tiled_f.m je ona koja reže na originalan ROI.
Funkcija 2: Cut_FixMPS_Tiled_f.m reže na smanjeni ROI koji je trebao riješiti problem,
ali na kraju ništa ne radi. Sve vrijednosti ispadaju NaN.
Ako hoces to isprobati samo sve funkcije 1 s funkcijama 2 u skripti 'Processing'.
