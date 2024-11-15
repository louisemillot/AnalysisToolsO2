#ifndef NCOLL_H
#define NCOLL_H

double Ncoll[100] = {
1983.83,    //  0-1%    
1868.86,    //  1-2%    
1770.18,    //  2-3%    
1679.54,    //  3-4%    
1596.57,    //  4-5%    
1520.04,    //  5-6%    
1448.43,    //  6-7%    
1381.50,    //  7-8%    
1318.23,    //  8-9%    
1258.54,    //  9-10%   
1201.89,    // 10-11%   
1147.87,    // 11-12%   
1096.98,    // 12-13%   
1048.45,    // 13-14%   
1001.63,    // 14-15%   
956.51,     // 15-16%   
913.45,     // 16-17%   
872.50,     // 17-18%   
833.08,     // 18-19%   
795.32,     // 19-20%   
759.33,     // 20-21%   
724.63,     // 21-22%   
691.21,     // 22-23%   
659.12,     // 23-24%   
628.15,     // 24-25%   
598.30,     // 25-26%   
569.58,     // 26-27%   
541.81,     // 27-28%   
515.15,     // 28-29%   
489.86,     // 29-30%   
465.72,     // 30-31%   
442.42,     // 31-32%   
419.79,     // 32-33%   
397.86,     // 33-34%   
376.91,     // 34-35%   
357.07,     // 35-36%   
338.06,     // 36-37%   
319.60,     // 37-38%   
301.77,     // 38-39%   
284.78,     // 39-40%   
268.61,     // 40-41%   
253.11,     // 41-42%   
238.26,     // 42-43%   
224.06,     // 43-44%   
210.46,     // 44-45%   
197.51,     // 45-46%   
185.21,     // 46-47%   
173.48,     // 47-48%   
162.24,     // 48-49%   
151.53,     // 49-50%   
141.39,     // 50-51%   
141.39,     // 50-51%   
131.83,     // 51-52%   
122.82,     // 52-53%   
114.25,     // 53-54%   
106.11,     // 54-55%   
98.39,      // 55-56%   
91.09,      // 56-57%   
84.23,      // 57-58%   
77.80,      // 58-59%   
71.74,      // 59-60%   
66.08,      // 60-61%   
60.82,      // 61-62%   
55.93,      // 62-63%   
51.35,      // 63-64%   
47.06,      // 64-65%   
43.03,      // 65-66%   
39.29,      // 66-67%   
35.85,      // 67-68%   
32.69,      // 68-69%   
29.76,      // 69-70%   
27.07,      // 70-71%   
24.59,      // 71-72%   
22.31,      // 72-73%   
20.23,      // 73-74%   
18.32,      // 74-75%   
16.59,      // 75-76%   
15.02,      // 76-77%   
13.57,      // 77-78%   
12.26,      // 78-79%   
11.06,      // 79-80%   
9.97,      // 80-81%   
8.99,      // 81-82%   
8.10,      // 82-83%   
7.29,      // 83-84%   
6.56,      // 84-85%   
5.90,      // 85-86%   
5.29,       // 86-87%   
4.75,       // 87-88%   
4.25,       // 88-89%   
3.81,       // 89-90%   
3.43,       // 90-91%   
3.09,       // 91-92%   
2.79,       // 92-93%   
2.54,       // 93-94%   
2.32,       // 94-95%   
2.13,       // 95-96%   
1.97,       // 96-97%   
1.83,       // 97-98%   
1.71,       // 98-99%   
1.59,       // 99-100%  
}

double GetMeanNcoll(float* centRange, TFile NcollFile) {
  int centLow = (int)centRange[0];
  int centHigh = (int)centRange[1];

  // calculate mean Ncoll in Ncoll[] table between centLow and centHigh; I think it's correct because centrality distribution of collisions is flat?
  double meanNcoll = 0;
  for(int iCent = centLow; iCent < centHigh; iCent++){
    meanNcoll += Ncoll[iCent];
  }
  meanNcoll = meanNcoll / (centHigh - centLow);

  return meanNcoll;
}


#endif