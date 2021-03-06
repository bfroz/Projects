# Artificial Crawlers

The current status: **Finished (published)**

### Abstract
Lung cancer is the major cause of death among patients with cancer throughout the world. The main symptom that indicate the lung cancer is the presence of lung nodules. This work proposes a methodology to classify lung nodule and non-nodule using texture features. The state-of-art of the presented work are the adaption of the Artificial Crawlers and Rose Diagram techniques for representing patterns over 3D images. 

Several information are extracted based on the texture behavior of these methods, allowing the correct classification of lung nodules candidates using Support Vector. Objective: This work proposes a methodology to classify lung nodule candidates and non-nodule candidates based on computed tomography (CT) images. Methodology: The Lung Image Database Consortium (LIDC-IDRI) image database is employed for our tests. 

Three techniques are employed to extract texture measurements. The first technique is artificial crawlers (ACs), an artificial life algorithm. The second technique uses the rose diagram (RD) to extract directional measurements. The third technique is a hybrid model that combines texture measurements from artificial crawlers and the rose diagram. The support vector machine (SVM) classifier with a radial basis kernel is employed. Results: In the testing stage, we used 833 scans from the LIDC-IDRI database. 

### Results

For the application of the methodology, we decided to divide the whole database into two groups, training and testing. We used partitions of training and testing of 20/80%, 40/60%, 60/40% and 80/20%. The division was repeated 5 times at random. We reached a mean accuracy (mACC) of 94.30%, a mean sensitivity (mSEN) of 91.86%, a mean specificity (mSPC) of 94.78%, a coefficient of accuracy variance (CAv) of 1.61% and a mean area under the receiver operating characteristic (mROC) curves of 0.922. 

#### Conclusion

Lung cancer has the highest mortality rate and one of the smallest survival rates after diagnosis. An early diagnosis increases the survival chance of patients. The proposed methodology is a useful tool for specialists in the detection of nodules. We believe we contribute for the expert system field because 1) the adaption of the Artificial Crawlers and Rose Diagram methods as 3D texture descriptors is innovative and contains great potential; 2) we adapted and developed measurements from the 3D texture descriptors; and 3) the simplicity and discriminative power of the methodology can be extended to applications based on images with other contexts.


For further details, please read the paper "Lung nodule classification using artificial crawlers, directional texture and support vector machine" or access the published paper at: https://doi.org/10.1016/j.eswa.2016.10.039.