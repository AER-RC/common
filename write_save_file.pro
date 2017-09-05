pro write_save_file, binTAPE, file_type=file_type

; binTAPE -- string, path to binary TAPE[11,12,13] file
; fType -- int, file type (radiance, transmittance, etc.; see doc in
;    /project/rc/rc2/mshep/idl/patbrown/read_lbl_file_dbl.pro)

if not keyword_set(file_type) then file_type=0

read_lbl_file, binTAPE, spectrum, wavenum, file_type=file_type, /double

; temporary IDL save file that can be read in and plotted with Python
save, filename='LBLRTM_output.sav', spectrum, wavenum

end
