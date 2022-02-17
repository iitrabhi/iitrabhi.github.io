---
layout: post
title: "Automated RCC design drawing checking system."
description: "Automating the RCC drawing checking procedure in VBA."
langs: [VBA]
year: "March 2015"
date: 2015-03-01
typora-root-url: ../../../../website
---

In the second year at NPCIL, I was finally transferred to the structural engineering department. One of the best professional experiences at NPCIL was in my first meeting with the consultants. The meeting was on the finalization of the depth of the raft foundation. The model was built in STAAD Pro, and the discussion was on the unrealistic bending moment in the model. 

During my M.Tech, I developed a particular liking to plates and shells because of our professor, Dr. Anupam Chakrabarti. He made the class fun, and I learned a lot from the course. During that time, I learned about moment-peaking, a side-effect of modeling columns as line elements and attaching them to plate elements at a single node. When I presented this information to the meeting attendees, they had a long discussion, and finally, everyone agreed with my point. This was my first meeting, and I had a win in my pocket.

Moreover, after the meeting the head of the Nuclear Division of L&T came to me and thanked me. He presented a calculation of the amount of money saved from this small decision. The total thickness of the raft foundation was reduced by 400mm, and it saved them almost 1.2 crore rupees 🤯. This experience developed a particular interest in the finite element method in me and led me to research more in this direction during my days at NPCIL.

It was also in the same meeting that I came to know about the working philosophy of NPCIL with consultants. We would go on such meetings to their office and check the engineering drawings developed by them. This was a tedious process and required manual checking of the design data in the reports and the design drawings. Again, I resorted to further enhancing my skills in VBA and developing a system for automated testing of the drawing and reports in Excel itself. 

Here are the highlights of that work.

- A system was developed in VBA-Excel to generate a database of all the sections used in the building. The user must enter the section name and its information in the database sheet. 
- Once the report data is entered into the database, we could pull the output forces and moments from STAAD Pro. This step was also automated in Excel with the help of VBA. STADD Pro allows us to connect its output database to Excel with the help of VBA. I used this connection to pull all of the data into Excel automatically.
- And finally, the last part was just to hit run and enjoy. All checks were carried out automatically, and a report was generated.

My boss appreciated this development, and he also suggested my name for the Young Scientist Award.
