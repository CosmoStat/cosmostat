ve description du package :
Les fonctions principales sont mrs_msvsts_IUWT_denoising (qui 
effectue le denoising sur une carte HEALPix à l'aide de la méthode MS-
VSTS + IUWT, avec des keywords pour faire de l'inpainting (mask) et 
de l'extraction du background (background)), et 
mrs_msvsts_curvelets_denoising (denoising MS-VSTS + curvelets).

Les autres fonctions sont des subroutines : 
mrs_msvsts_IUWT_param_computing calculee les paramètres de la VST en 
racine carrée pour chaque échelle de la IUWT. Ce calcul est utilisé 
dans toutes les autres routines.
mrs_msvsts_IUWT_transform effectue la transformation MS-VSTS+IUWT sur 
une carte HEALPix. En en sortie, on a la transformée en ondelettes 
"redressée" (structure IDL semblable à la sortie de mrs_wttrans).
mrs_msvsts_IUWT_hypothesis_testing effectue la transformation MS-VSTS
+IUWT sur une carte HEALPix en appelant mrs_msvsts_IUWT_transform, 
puis effectue des tests d'hypothèses pour déterminer les coefficients 
significatifs. La sortie est le support de multi-résolution.
Une subroutine mrs_msvsts_IUWT_reconstruction (dans le fichier 
mrs_msvsts_IUWT_denoising) effectue la reconstruction itérative de 
l'estimée finale, à partir des données et du support multi-résolution.


