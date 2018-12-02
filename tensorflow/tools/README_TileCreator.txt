
TileCreator
Tool to create tiles from 2D or 3D simulation data (or other grid data) with different channels.
Includes data augmentation: scaling, shifting and rotation. scaling and rotation is also applied to vector data in channels.

optional second set of sturctured data with matching dimensionality and size compatible to the main data

arbitrary label alongside the structured data


0. Dependencies
python v3.6
numpy v1.14.0 (v1.15.4 for scaled data)
scipy v0.18.1
imageio v2.1.2

1. Setup
import tilecreator_t as tc
TC = tc.TileCreator(tileSize, simSize=64, dim=2, densityMinimum=0.02, channelLayout_main=C_LAYOUT['dens_vel'], useScaledData=True, channelLayout_scaled=C_LAYOUT['dens'], scaleFactor=2, useDataBlocks=False, useLabels=False, partTrain=0.8, partTest=0.2, partVal=0, logLevel=LOG_LEVEL_WARNING)

tileSize: size of the tiles to create. must be less or equal the simulation size. Assumed to have enough bounds for rotation augmentation if active. (previously lowTileSize)
simSize: size of the input simulation data. (previously lowSimSize)
dim: dimension of main (and scaled) data. can be 2 or 3.
densityMinimum: minimum avg. density in a tile. To prevent generating empty tiles
channelLayout_main: what type of data the different channels contain as a comma separaed sring of channel keys. used in augmentation. Keys:
	d: default. not augmented
	v: vector data. needs 2 (2D; x,y) or 3 (2D, 3D; x,y,z) components with matching labels. the order does not matter and one empty label is accepted. format: v[label](x|y|z). examples: '...,vx,vy,vz,...', 'd,vVELx,vVELy,d,vVortz,vVortx,d,vVorty'

useScaledData: an optional second structured dataset with a fixed scaling factor to the main data. will be augmented to fit the augmentation of the main data (e.g. same rotation)
channelLayout_scaled: same as channelLayout_main for active scaled data
scaleFactor: the scaling factor between main and scaled data. must be larger or equal to 1. must be a multiple of 1/16 to avoid numerical errors. simSize*scaleFactor and tileSize*scaleFactor must be whole numbers.

useDataBlocks: an optional grouping of data using block ids. the ids have to be provided when adding data to the tilecrator. enables the creation of blocks of tiles with matching augmentation (i.e. like an additional dimension that is not agmented). can be used to create augmented (time-) sequences of data
useLabels: an optional set of data that is not augmented. no type or structure is assumed
partTrain, partTest (, partVal): the relative sizes of data sets for train and testing mode (machine learing). val data is currently unused and inaccessible and should be left at 0.
logLevel: how much information to print. tc.LOG_LEVEL_ERROR, tc.LOG_LEVEL_WARNING, tc.LOG_LEVEL_INFO, tc.LOG_LEVEL_DEBUG.

1.1 Setup Data Augmentation
TC.initDataAugmentation(rot=2, minScale=0.85, maxScale=1.15 ,flip=True)
Enables data augmentation for main (and scaled) data with the specified parameters. Set augment=True when calling TC.selectRandomTiles() to create augmented data.

rot: type of rotation augmentation (NOT an angle limit). 1: fixed 90° roatations, 2: full rotation (assumes enough space for boundaries in data), else: no rotation)
minScale, maxScale: limits for scaling. set both to 1 to disable scaling.
flip: flipping (mirroring) of data

tile size is checked to be compatible with the specified rotation and scaling parameters.

1.2 Using block data from fluiddataloader
The fluiddataloader merges blocks into channels. Use tc.blockFromChannelsToSequence(data, blockLength) to convert the data to the required format.

2. Adding data
TC.addData(main, scaled=None, labels=None, blocks=None)
Adds data to the TileCreator. The data is split into training,validation and testing set afterwards. Data can be cleared with TC.clearData()

main: the main data, must match simSize and channels specified in the constructor
scaled: required when useScaledData==True, ignored otherwise. must match simSize*scaleFactor and channels specified in the constructor. must be the same amount as the main data.
labels: required when useLabels==True, ignored otherwise. list/iterable of arbitrary data. must be the same amount as the main data.
blocks: required when useDataBlocks==True, ignored otherwise. list of block ids (int). can be unsorted. will be sorted according to the id, the order of data with the same id (within the same block) is preserved.

3. Batch creation
TC.selectRandomTiles(selectionSize, isTraining=True, augment=False, blockSize = 1, squeezeZ=False, squeezeBlocks=True)

selectionSize: number of tiles to create
isTraining: whether to use data from the training or testing set
augment: whether to augment the data. requires data augmentation to be initialized (TC.initDataAugmentation).
blockSize: what block size to use if block data is active. ignored otherwise
squeezeZ: whether to squeeze/collapse the z dimension/axis of main (and scaled) data when using 2D data.
squeezeBlocks: whether to squeeze/collapse the block dimension/axis of main (and scaled) when block size is 1 or block data is inactive

returns:
main[,scaled][,labels]: main and scaled are np.ndarray with shape: n[,b][,z],y,x,c with z,y,x mathching the (scaled) tile size and c channels. labels is a list

4. Image output
tc.savePngs(tiles, path, tile_format='NYXC', imageCounter=0, tiles_in_image=[1,1], channels=[0], save_gif=False, plot_vel_x_y=False, save_rgb=None, rgb_interval=[-1,1])

tiles: tile data to save as images (.png).
Path: directory to write images to.
tile_format: data format of the imput tile. valid formats are 'yx','yxc','nyx','nyxc','byxc', 'nbyxc', not case sensitive.
imageCounter: as this method can output multiple images for sequences of tiles, the images are numbered starting with this number.
tiles_in_image: number of tiles along x and y dimension to combine to a single image.
channels: channels to save as grayscale image.
save_gif:
plot_vel_x_y: save 2nd and 3rd channel as vector quiver plot.
save_rgb: list of lists of length 2 or 3. save these channels as RGB image.
rgb_interval: data values are mapped to [0,1] from this interval when writing RGB images.

returns: next image counter (imageCounter + number of images written)


5. known issues
very small tiles might cause errors when cutting tiles during augmentation due to rounding.