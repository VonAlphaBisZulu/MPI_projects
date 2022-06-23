Zur Verwendung wird eine angepasste Version der Datei "MCS_MILP" benötigt, da der Wert obj.b_eq in manchen Berechnungen
angepasst wird.
%%%%%Verbesserung der Big-M-Formulierung über einen iterative Ansatz %%%%%

Das Skript "AMIP_ALP" führt den Algorithmus. Es stehen die Netze Tiny Net oder Small Example 3 zur Verfügung.
Zum Ändern des Netzes muss das gewünschte jeweils auskommentiert werden. 

Die Netzwerk-Dateien in dem Ordner sind bereits in der Form, wie sie nach der ersten Berechnung (FVA-Ansatz) in computeM vorliegen.
Solver: Intlinprog
%%%%% Iterative Domain Reduction %%%%

Ausführung des Algorithmus durch "iterative_domain_reduction". Es stehen die Netze Tiny Net oder Small Example 3 zur Verfügung.
Zum Ändern des Netzes muss das gewünschte jeweils auskommentiert werden.  Für jedes Netz wurden ebenfalls zwei Startwerte für M genutzt. 
Entsprechend muss der Wert kommentiert/auskommentiert werden. 

Die Netzwerk-Dateien in dem Ordner sind bereits in der Form, wie sie nach der ersten Berechnung (FVA-Ansatz) in computeM vorliegen.
Solver: Intlinprog

%%%% Abschätzung Big-M-Werte über Skalierung %%%
Bei diesem Ansatz wird b_f  der Farkas-Nebenbedinung geändert. Für jede Berechnung existiert ein Ordner mit dem jeweiligem Skript, der computeM 
Funktion und dem Netz. Die Berechnungen für CPLEX und Intlinprog werden direkt in einem Skript ausgeführt.

%%% Optimale Big-M-Werte %%%

Im Ordner "Berechnung MCS Eckenenumeration" werden die MCS für das Small Example 3 und das Tiny-Net durch Bestimmung der M_i^{opt} mittels Ecken-
enumeration (MPT3) durchgeführt. 

https://www.mpt3.org/Main/Installation
Solver: CPLEX

Im Ordner "Berechnung MCS Eckenenumeration2" werden die MCS für das Small Example 3 und das Tiny-Net durch Bestimmung der M_i^{opt} mittels Ecken-
enumeration (MPT3) durchgeführt. Zusätzlich werden die berechneten M_i{opt} um den Faktor 0.9 verringert.

https://www.mpt3.org/Main/Installation
Solver: CPLEX

Im Ordner "Vergleich M_iopt und Ecken" führt das Skript "script_vertices_milp" zunächst eine Berechnung der M_i^{opt} mittels
Eckenenumeration aus und anschließend werden die M_i^{opt} erneut über ein MILP mittels Indicator-Constraints bestimmt. 
Es stehen das Tiny-Net und Small Example 3 zur Verfügung und können durch auskommentieren gewählt werden.

Solver: CPLEX

Der Ordner "Berechnung ECC2Comp" beinhaltet das Skript "Skript_ECC2Comp". Damit wird eine MCS Berechnung gestartet, wobei das MILP für die M_i^{opt}
bis zu einem Knotenlimit von 50k berechnet werde. Die Überprüfung für \Tilde(b) ist zunächst, aufgrund zu langer Laufzeit ausgeklammert. Im Allgemeinen erfüllt
ein randomisierter Vektor die Ansprüche.  

Solver: CPLEX

Der Ordner "Berechnung ECC2Comp2" fügt den approximierten M_i^{opt} einen Offset hinzu. 

Solver CPLEX

