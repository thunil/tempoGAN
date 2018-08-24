import os
#os.environ["CUDA_VISIBLE_DEVICES"]="0"

#using density as input of generator, full tempoD
#2D
os.system('python tempoGAN.py randSeed 42 out 1 trainingIters 40000 lambda 5.0 lambda2 -0.00001 discRuns 2 genRuns 2 alwaysSave 1 fromSim 2006 toSim 2006 outputInterval 200 genValiImg 1 dataDim 2 batchSize 16 useVelocities 0 useVorticities 0 gif 0 genModel gen_resnet discModel disc_binclass basePath ../2ddata_gan/ loadPath ../2ddata_sim/ lambda_t 0.0 lambda_t_l2 0.0 frameMax 120 data_fraction 1.0 adv_flag 0 dataAugmentation 0 rot 2 decayLR 1 load_model_test 0 load_model_no 199 simSize 128 tileSize 128')
#3D
os.system('python tempoGAN.py randSeed 42 out 1 trainingIters 40000 lambda 5.0 lambda2 -0.00001 discRuns 2 genRuns 2 alwaysSave 1 fromSim 3006 toSim 3006 outputInterval 200 genValiImg 1 dataDim 3 batchSize 1 useVelocities 0 useVorticities 0 gif 0 genModel gen_resnet discModel disc_binclass basePath ../3ddata_gan/ loadPath ../3ddata_sim/ lambda_t 0.0 lambda_t_l2 0.0 frame_min 100 frameMax 110 data_fraction 1.0 adv_flag 0 dataAugmentation 0 rot 2 decayLR 1 load_model_test 0 load_model_no 34 simSize 64 tileSize 35 overlap 5')


#using density and velocity as inputs of generator, full tempoD
os.system('python tempoGAN.py randSeed 42 out 1 trainingIters 40000 lambda 5.0 lambda2 -0.00001 discRuns 2 genRuns 2 alwaysSave 1 fromSim 2006 toSim 2006 outputInterval 200 genValiImg 1 dataDim 2 batchSize 16 useVelocities 1 useVorticities 0 gif 0 genModel gen_resnet discModel disc_binclass basePath ../2ddata_gan/ loadPath ../2ddata_sim/ lambda_t 0.0 lambda_t_l2 0.0 frameMax 120 data_fraction 1.0 adv_flag 0 dataAugmentation 0 rot 2 decayLR 1 load_model_test 0 load_model_no 199 simSize 128 tileSize 128')

#using density, velocity and vorticity as inputs of generator, full tempoD
os.system('python tempoGAN.py randSeed 42 out 1 trainingIters 40000 lambda 5.0 lambda2 -0.00001 discRuns 2 genRuns 2 alwaysSave 1 fromSim 2006 toSim 2006 outputInterval 200 genValiImg 1 dataDim 2 batchSize 16 useVelocities 1 useVorticities 1 gif 0 genModel gen_resnet discModel disc_binclass basePath ../2ddata_gan/ loadPath ../2ddata_sim/ lambda_t 0.0 lambda_t_l2 0.0 frameMax 120 data_fraction 1.0 adv_flag 0 dataAugmentation 0 rot 2 decayLR 1 load_model_test 0 load_model_no 199 simSize 128 tileSize 128')

#using density, velocity and vorticity as inputs of generator, l2 loss for temporal
os.system('python tempoGAN.py randSeed 42 out 1 trainingIters 40000 lambda 5.0 lambda2 -0.00001 discRuns 2 genRuns 2 alwaysSave 1 fromSim 2006 toSim 2006 outputInterval 200 genValiImg 1 dataDim 2 batchSize 16 useVelocities 1 useVorticities 1 gif 0 genModel gen_resnet discModel disc_binclass basePath ../2ddata_gan/ loadPath ../2ddata_sim/ lambda_t 0.0 lambda_t_l2 1.0 frameMax 120 data_fraction 1.0 adv_flag 0 dataAugmentation 0 rot 2 decayLR 1 load_model_test 0 load_model_no 199 simSize 128 tileSize 128')
