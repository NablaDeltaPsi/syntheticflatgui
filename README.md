
## SyntheticFlatGUI <img src="SyntheticFlatGUI.ico" width=30 align="right">

Berechnen eines synthetischen Flats zur Vignettierungs-Korrektur 

<img src="readme_images/GUI.png">

### Start

**Mit Python**  
Für alle im Header von **SyntheticFlatGUI.py** aufgeführten, benötigten Pakete, Prüfen ob sie in der Python Distribution installiert sind (`python -m pip show <package>`) und Installieren wenn nicht (`python -m pip install <package>`). Wenn alle Voraussetzungen erfüllt sind, mit `python SyntheticFlatGUI.py` ausführen.

**Executable für Windows-Benutzer**  
Für Windows-Benutzer ist unter Releases eine ZIP-Datei mit kompiliertem Programm verfügbar. Herunterladen, Entzippen und **SyntheticFlatGUI.exe** ausführen.

### Einleitung
Um in der Astrofotografie die Vignettierung von Objektiven zu Korrigieren, muss vor dem Stacken ein (Master-) Flat-Frame von den Einzelbildern dividiert werden. Normalerweise erstellt man das Master-Flat, indem mehrere Bilder einer weißen Wand aufgenommen und kombiniert werden, um möglichst wenig zusätzliches Rauschen ins Bild zu bringen.

Dennoch bringen qualitativ minderwertige oder unpassende Flats oft auch Probleme wie Raining Noise, Überkorrektur oder ringartige Strukturen ins Bild. Dieses Programm ermöglicht die einfache Erstellung synthetischer Flats, um diese Probleme zu umgehen.

### Vorgehen
Die synthetischen Flats können entweder aus üblichen manuellen Wand-Flats oder sogar aus den Lights berechnet werden, solange keine größeren flächigeren Strukturen zu sehen sind (Milchstraße). Zur Erstellung wird das Bild erst vom Bias/Offset-Wert befreit, der manuell angegeben oder aus einem dunklen Bias-Frame geladen werden kann. Anschließend wird ein radiales Profil berechnet, mit einer auswählbaren Sigma-Clipping Statistik manipuliert um Spitzen von Sternen zu entfernen, und anschließend weiter geglättet (Savitzki-Golay Filter). Ausgegeben wird ein TIF in der gleichen Größe wie das Original.

<img src="readme_images/combined.png">

### Verwendung
Zur Verwendung sollte von den Lights ebenfalls erst der gleiche Bias-Wert abgezogen werden. Anschließend kann das Ergebnis durch das synthetische Flat dividiert werden, um ein korrigiertes Bild zu erhalten. In Siril sind bspw. die entscheidenden Kommandos zur Verarbeitung des Flats:

```
# tif to fits
load masters/master_flat.tif
save masters/master_flat
```

```
# convert with synthetic bias and hotpixel list
calibrate light -bias="=500" -flat=../masters/master_flat -cc=bpm ../masters/master_bias_hotpixels.lst -cfa -equalize_cfa -debayer -prefix=cal_
```


### Schaltflächen
- Load files  
Auswählen eines oder mehrerer RAW-Bilder.
- Set bias value  
Abziehen eines Kameraabhängigen Offset-Wertes vor Berechnung aller folgenden Funktionen
- Bias from file  
Auswählen eines dunklen Bias Frames zur Berechnung des Bias-Wertes. Zur Berechnung wird ein 2-Sigma-Clipping verwendet, cold- und hotpixel werden also ignoriert.
- Start  
Anwenden der Funktionen unter "Options" 

### Options
Ein/Aus Schalter der Hauptfunktionen
- Correct gradient  
Ein Gradient im Bild muss vor der Erstellung eines Flats korrigiert werden. Ansonsten wird das Sigma-Clipping aushebelt und es kommt zu Unregelmäßigkeiten im radialen Profil!
- Calculate histogram  
Speichert eine CSV mit dem Histogramm des Rohbilds (RGB) ab.
- Calculate radial profile  
Ermittelt die Helligkeitskurve des Bildes in Abhängigkeit zum Abstand zur Bildmitte und speichert CSVs dazu ab. Zur Ermittlungslogik siehe "Statistics".
- Export synthetic flat  
Mit dem berechneten radialen Profil (Option wird automatisch gechecked) wird ein normiertes Flat (16-bit TIF) berechnet, das gleich groß ist wie das Ursprungsbild.

### Settings
Kleinere Schalter zum Beeinflussen der Hauptfunktionen
- Write pickle file  
Schreiben einer PKL Datei des gedebayerten Bildes zum schnelleren Ausführen darauffolgender Läufe
- Histogram of largest circle  
Das Histogramm eines Flats wird normalerweise durch das rechteckige Beschneiden beeinflusst. Mit dieser Option wird das Histogramm nur für den größtmöglichen Kreis im Bild berechnet.
- Extrapolate inside max  
Radiale Profile können ihr Maximum statt bei Radius 0 bei größeren Radii aufweisen, wodurch im synthetischen Flat Ringe entstehen würden. Mit dieser Option wird das radiale Profil beim Maximum abgeschnitten und zum Zentrum hin mit einer quadratischen Funktion extrapoliert.
- Export corrected input images  
Es wird nicht nur das synthetische Flat als TIF ausgegeben, sondern auch das Originalbild, das gradientenbereinigte Originalbild sowie ein flat-bereinigtes Originalbild.
- Grey synthetic flat  
Für das synthetische Flat wird für alle vier RGGB Untergitter das gleiche (grüne) radiale Profil verwendet. Im Ergebnisbild ist dann kein Karomuster zu sehen.

### Statistics
Auswählen der verwendeten Statistik zur Berechnung des radialen Profils. Auf einem Ring (vgl. Bild oben) liegen mehrere Pixel mit unterschiedlichen Werten. Aus diesen Werten pro Ring kann der Mittelwert, der Median, das Maximum, oder das Minimum genommen werden, oder sie werden erst mit einem Sigma-Clipping aussortiert (empfohlen) und erst dann der Mittelwert genommen.

### Resolution
Die Berechnung des radialen Profils für ein 24 MP Bild kann einige Zeit dauern, obwohl eine so große Auflösung gar nicht nötig ist. Deshalb kann mit dieser Option die Statistik mit einem verkleinerten Bild berechnet werden. Wichtig: Die Verkleinerung spielt nur für die statistischen Funktionen (Pixelmap, Histogramm, radiales Profil) eine Rolle, alle ausgegebenen Bilder haben aber immer stets die gleiche Größe wie das Ursprungsbild! 

### Troubleshooting
- Bisher nur mit Sony ARW Bildern getestet. Bei Problemen mit anderen Formaten mir gerne Testbilder schicken.
- Bisher nur auf RGGB-Pattern ausgelegt. Bei Bedarf baue ich gerne eine Option für andere Pattern ein.
- Mit zu hohem Bias-Wert kann es zu Problemen kommen. Eventuell mal ohne Bias (0) testen. Als Richtwert: Bei meiner Sony A7 II liegt der Bias bei 512.
- Ein Gradient im Bild hebelt das Sigma-Clipping der radialen Profile aus. Bei unschönen Profilen unbedingt "Correct gradient" einschalten!

