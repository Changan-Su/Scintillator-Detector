/public1/home/a8s001349/soft/geant4-11.2.2/install/share/Geant4/data

DATA_ROOT="/public1/home/a8s001349/soft/geant4-11.2.2/install/share/Geant4/data"

# 小工具：按前缀匹配一个子目录（取第一个匹配）
finddir() { ls -d "${DATA_ROOT}/$1"* 2>/dev/null | head -n 1; }

export G4ENSDFSTATEDATA=$(finddir G4ENSDFSTATE)         # 例如 .../G4ENSDFSTATE2.3
export G4LEDATA=$(finddir G4EMLOW)                      # 例如 .../G4EMLOW8.5
export G4LEVELGAMMADATA=$(finddir PhotonEvaporation)    # 例如 .../PhotonEvaporation5.7
export G4RADIOACTIVEDATA=$(finddir RadioactiveDecay)    # 例如 .../RadioactiveDecay5.6
export G4NEUTRONHPDATA=$(finddir G4NDL)                 # 例如 .../G4NDL4.7
export G4ABLADATA=$(finddir G4ABLA)                     # 例如 .../G4ABLA3.1
export G4PARTICLEXSDATA=$(finddir G4PARTICLEXS)     # 例如 .../G4PARTICLEXSDATA3.1
export G4PIIDATA=$(finddir G4PII)                       # 例如 .../G4PII1.3

# 打印确认
echo "== Geant4 DATA =="
env | grep -E 'G4.*DATA' | sort
