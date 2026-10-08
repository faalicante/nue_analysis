set basePath to "/Users/fabioali/SND@LHC/nue_analysis/CNN/output/presentazione_progetto/"
tell application "Keynote"
  set projectDeck to open POSIX file (basePath & "Sintesi_progetto_CNN_SND.pptx")
  save projectDeck in POSIX file (basePath & "Sintesi_progetto_CNN_SND.key")
  export projectDeck to POSIX file (basePath & "Sintesi_progetto_CNN_SND.pdf") as PDF
  return count of slides of projectDeck
end tell
