Pokretanje skripti:
1. Capturing
2. Unpacking
3. Processing

Ono na �ta treba� paziti je raspored foldera.
Nisam nigdje napravio to da automatski radi foldere, pa se to mora ru�no podesiti. 
U tree.txt imas prikazano cijelo stablo. Znam da bi bilo lak�e da je to automatski,
ali sad nemam vremena za to. Jo� samo napomena, to stablo ide od .\data\Scenes\

Kao �to je najavljeno, MPS slike se spremaju u folder, ostale idu u �elije ili matrice.
Rezultati se spremaju isto u �eliju. 
Nisam probao skriptu 'Processing' za sve slike, samo za jednu. Sutra �u to izvrtiti.

Devignetting jo� nije popravljen. 
Funkcija 1: Cut_Tiled_f.m je ona koja re�e na originalan ROI.
Funkcija 2: Cut_FixMPS_Tiled_f.m re�e na smanjeni ROI koji je trebao rije�iti problem,
ali na kraju ni�ta ne radi. Sve vrijednosti ispadaju NaN.
Ako hoces to isprobati samo sve funkcije 1 s funkcijama 2 u skripti 'Processing'.
