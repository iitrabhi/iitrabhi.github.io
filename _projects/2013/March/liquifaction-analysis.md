---
layout: post
title: "Automated liquifaction analysis"
description: "Creating a Excel app to automate the procedure of liquifaction analysis."
langs: [VBA]
year: "March 2013"
typora-root-url: ../../../../website
---

Even though I had a degree in structural engineering, somehow, I got placed in the Geotech department instead of the structural at my job. The first project I was given was to carry out the liquefaction analysis for one of the Nuclear power plant sites. The current method was to fill out data in hundreds of sheets in excel and then calculate the results from each sheet individually; this was such an unoptimized way of doing things. It took almost a month to enter the data and carry out the complete analysis, and then prepare a report out of it.

This was an opportunity to showcase my coding competence and learn a new language - VBA. So I did that. Here are the highlights of that work.

- All of the approximate 150 sheets were converted into a three-sheet model. The first sheet consisted of a database of all the 150 sheet data. The second sheet was the calculation engine, and the third was the output report.
- The logic of connecting the data to the calculation engine to the output report was done in VBA. I learned VBA from an excellent book by [Shepherd Richard](https://www.flipkart.com/excel-2007-vba-macro-programming/p/itmetj5vdzfuczbj).
- It took me about a month to learn and develop the whole system. But, once it was built, the analysis time was reduced from one month to around 5min for each run. And as the database was centralized, there was no chance of human error.
- I was proud of this work.
