#####batch process

#0.05

python3 main.py 40 8 /home/ch/Documents/2_FlowPyTesting/04/005/cl1 /home/ch/Documents/qgis/PROCESSING/01_DATA/DEM10m_FINSTERTAL_resized.tif /home/ch/Documents/qgis/PROCESSING/01_DATA/005/cl1_PRA005.tif flux=0.003 max_z=270

python3 main.py 35 8 /home/ch/Documents/2_FlowPyTesting/04/005/cl2 /home/ch/Documents/qgis/PROCESSING/01_DATA/DEM10m_FINSTERTAL_resized.tif /home/ch/Documents/qgis/PROCESSING/01_DATA/005/cl2_PRA005.tif flux=0.003 max_z=270

python3 main.py 30 8 /home/ch/Documents/2_FlowPyTesting/04/005/cl3 /home/ch/Documents/qgis/PROCESSING/01_DATA/DEM10m_FINSTERTAL_resized.tif /home/ch/Documents/qgis/PROCESSING/01_DATA/005/cl3_PRA005.tif flux=0.003 max_z=270

#0.25


python3 main.py 40 8 /home/ch/Documents/2_FlowPyTesting/04/025/cl1 /home/ch/Documents/qgis/PROCESSING/01_DATA/DEM10m_FINSTERTAL_resized.tif /home/ch/Documents/qgis/PROCESSING/01_DATA/025/cl1_PRA025.tif flux=0.003 max_z=270

python3 main.py 35 8 /home/ch/Documents/2_FlowPyTesting/04/025/cl2 /home/ch/Documents/qgis/PROCESSING/01_DATA/DEM10m_FINSTERTAL_resized.tif /home/ch/Documents/qgis/PROCESSING/01_DATA/025/cl2_PRA025.tif flux=0.003 max_z=270

python3 main.py 30 8 /home/ch/Documents/2_FlowPyTesting/04/025/cl3 /home/ch/Documents/qgis/PROCESSING/01_DATA/DEM10m_FINSTERTAL_resized.tif /home/ch/Documents/qgis/PROCESSING/01_DATA/025/cl3_PRA025.tif flux=0.003 max_z=270

#0.45

python3 main.py 40 8 /home/ch/Documents/2_FlowPyTesting/04/045/cl1 /home/ch/Documents/qgis/PROCESSING/01_DATA/DEM10m_FINSTERTAL_resized.tif /home/ch/Documents/qgis/PROCESSING/01_DATA/045/cl1_PRA045.tif flux=0.003 max_z=270

python3 main.py 35 8 /home/ch/Documents/2_FlowPyTesting/04/045/cl2 /home/ch/Documents/qgis/PROCESSING/01_DATA/DEM10m_FINSTERTAL_resized.tif /home/ch/Documents/qgis/PROCESSING/01_DATA/045/cl2_PRA045.tif flux=0.003 max_z=270

python3 main.py 30 8 /home/ch/Documents/2_FlowPyTesting/04/045/cl3 /home/ch/Documents/qgis/PROCESSING/01_DATA/DEM10m_FINSTERTAL_resized.tif /home/ch/Documents/qgis/PROCESSING/01_DATA/045/cl3_PRA045.tif flux=0.003 max_z=270

#0.65

python3 main.py 40 8 /home/ch/Documents/2_FlowPyTesting/04/065/cl1 /home/ch/Documents/qgis/PROCESSING/01_DATA/DEM10m_FINSTERTAL_resized.tif /home/ch/Documents/qgis/PROCESSING/01_DATA/065/cl1_PRA065.tif flux=0.003 max_z=270

python3 main.py 35 8 /home/ch/Documents/2_FlowPyTesting/04/065/cl2 /home/ch/Documents/qgis/PROCESSING/01_DATA/DEM10m_FINSTERTAL_resized.tif /home/ch/Documents/qgis/PROCESSING/01_DATA/065/cl2_PRA065.tif flux=0.003 max_z=270

python3 main.py 30 8 /home/ch/Documents/2_FlowPyTesting/04/065/cl3 /home/ch/Documents/qgis/PROCESSING/01_DATA/DEM10m_FINSTERTAL_resized.tif /home/ch/Documents/qgis/PROCESSING/01_DATA/065/cl3_PRA065.tif flux=0.003 max_z=270

print('HELL YEAH')

###calculate 0.25 cl3 with 35 and 40

python3 main.py 35 8 /home/ch/Documents/2_FlowPyTesting/04/025/cl3_35 /home/ch/Documents/qgis/PROCESSING/01_DATA/DEM10m_FINSTERTAL_resized.tif /home/ch/Documents/qgis/PROCESSING/01_DATA/025/cl3_PRA025.tif flux=0.003 max_z=270

python3 main.py 40 8 /home/ch/Documents/2_FlowPyTesting/04/025/cl3_40 /home/ch/Documents/qgis/PROCESSING/01_DATA/DEM10m_FINSTERTAL_resized.tif /home/ch/Documents/qgis/PROCESSING/01_DATA/025/cl3_PRA025.tif flux=0.003 max_z=270

python3 print('HELL YEAH')

###pra FOREST vs noFOREST

python3 main.py 30 8 /home/ch/Documents/2_FlowPyTesting/05/FOREST /home/ch/Documents/qgis/PROCESSING/01_DATA/FORESTAREA/DEM10m_FOREST.tif /home/ch/Documents/qgis/PROCESSING/01_DATA/FORESTAREA/PRA_FOREST.tif flux=0.003 max_z=270

python3 main.py 30 8 /home/ch/Documents/2_FlowPyTesting/05/noFOREST /home/ch/Documents/qgis/PROCESSING/01_DATA/FORESTAREA/DEM10m_FOREST.tif /home/ch/Documents/qgis/PROCESSING/01_DATA/FORESTAREA/PRA_noFOREST.tif flux=0.003 max_z=270

python3 print('HELL YEAH')


