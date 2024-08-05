# trackGrainsPF

Varun Srinivas, Zachary Croft, Marcel Chlupsa  
University of Michigan (2024)  

# Workflow

1. Load data.  
Use your own:
```Matlab
load('myData.mat')
```
Or the data provided here:
```
https://drive.google.com/file/d/1cMUoTJSyItMCOXVXclEm0irjinZgFZQA/view?usp=sharing
```
```Matlab
load('trackGrainsPF_data.mat')
```

2. Segment each 4D (3D space, 1D order parameters) map of order parameters.
```Matlab
[map_t250,numElement_t250] = segment_from_order_parameters(OP_t250,boundsmask);
[map_t900,numElement_t900] = segment_from_order_parameters(OP_t900,boundsmask);
```

3. Track grain id's between states
```Matlab
[tracked_grains] = track_grains_improved(map_t250,map_t900);
```

4. Take a look at the result.  
Open the tracked_grains variable and visualize a slice from each state:
```Matlab
figure, imshow(squeeze(map_t250(:,200,:)),[])
figure, imshow(squeeze(map_t900(:,200,:)),[])
```
![image](https://github.com/user-attachments/assets/5249a2c8-5eac-46bf-be8b-ab6394d92151)

THE END
