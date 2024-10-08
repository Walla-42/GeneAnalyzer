<!-- Improved compatibility of back to top link: See: https://github.com/othneildrew/Best-README-Template/pull/73 -->
<a name="readme-top"></a>
<!--
*** Thanks for checking out the Best-README-Template. If you have a suggestion
*** that would make this better, please fork the repo and create a pull request
*** or simply open an issue with the tag "enhancement".
*** Don't forget to give the project a star!
*** Thanks again! Now go create something AMAZING! :D
-->



<!-- PROJECT SHIELDS -->
<!--
*** I'm using markdown "reference style" links for readability.
*** Reference links are enclosed in brackets [ ] instead of parentheses ( ).
*** See the bottom of this document for the declaration of the reference variables
*** for contributors-url, forks-url, etc. This is an optional, concise syntax you may use.
*** https://www.markdownguide.org/basic-syntax/#reference-style-links
-->
[![Contributors][contributors-shield]][contributors-url]
[![Forks][forks-shield]][forks-url]
[![Stargazers][stars-shield]][stars-url]
[![Issues][issues-shield]][issues-url]
[![MIT License][license-shield]][license-url]
[![LinkedIn][linkedin-shield]][linkedin-url]



<!-- PROJECT LOGO -->
<br />
<div align="center">
  <a href="https://github.com/walla-42/GeneAnalyzer">
    <img src="Images/GeneAnalyzer.png" alt="Logo" width="1000" height="200">
  </a>

<h3 align="center>Gene Analysis Program</h3>

  <p align="center">
 I made this program when I first started learning python. It started as a way for me to become familiar with bioinformatic algorithms and loops in python and has been the inspiration for my Gene finder program and many other ideas I have had over the past two years. I am currently updating the program to utilize an OOP format so that it will be easier to add features in the future. Feel free to check out both programs to see what has changed!
    <br />
    <a href="https://github.com/walla-42/GeneAnalyzer"><strong>Explore the docs »</strong></a>
    <br />
    <br />
    <a href="https://github.com/walla-42/GeneAnalyzer">View Demo</a>
    ·
    <a href="https://github.com/walla-42/GeneAnalyzer/issues/new?labels=bug&template=bug-report---.md">Report Bug</a>
    ·
    <a href="https://github.com/walla-42/GeneAnalyzer/issues/new?labels=enhancement&template=feature-request---.md">Request Feature</a>
  </p>
</div>



<!-- TABLE OF CONTENTS -->
<details>
  <summary>Table of Contents</summary>
  <ol>
    <li>
      <a href="#about-the-project">About The Project</a>
      <ul>
        <li><a href="#built-with">Built With</a></li>
      </ul>
    </li>
    <li>
      <a href="#getting-started">Getting Started</a>
      <ul>
        <li><a href="#prerequisites">Prerequisites</a></li>
        <li><a href="#installation">Installation</a></li>
      </ul>
    </li>
    <li><a href="#usage">Usage</a></li> 
    <li><a href="#roadmap">Roadmap</a></li>
    <li><a href="#contributing">Contributing</a></li>
    <li><a href="#license">License</a></li>
    <li><a href="#contact">Contact</a></li>
    <li><a href="#acknowledgments">Acknowledgments</a></li>
  </ol>
</details>



<!-- ABOUT THE PROJECT -->
## About The Project

<!--[![Product Name Screen Shot][product-screenshot]]


<p align="right">(<a href="#readme-top">back to top</a>)</p> -->

### Built With

* [![Python][Python.js]][Python-url]


<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- GETTING STARTED -->
## Getting Started



### Prerequisites

1. Before you can run the program make sure you have an updated version of Python3 and BioPython installed on your system.


    ```sh
    pip install BioPython numpy matplotlib seaborn
    ```

### Installation

1. Clone the repository
   ```sh
   git clone https://github.com/walla-42/GeneAnalyzer.git
   ```

<p align="right">(<a href="#readme-top">back to top</a>)</p> 



<!-- USAGE EXAMPLES -->
<!--## Usage

To use this program you need to first clone the repository on your system and open the folder in your IDE. 
Navigate to the 'Main.py' file and hit 'Run'. The program will prompt you to enter a query. You can put in any protein target you wish. The program will search the Chembl database and return a list of identified targets based off of your search query. If you want to try some out, you can enter in queries like - "Coronavirus", "acetylcholinesterase", "LRRK2", or, "SV2C". Select a drug target from the list by typing the index of your desired protein target. 

"""NOTE: THE PROGRAM WILL NOT PERFORM WELL IF YOU CHOOSE A QUERY THAT IS NOT IDENTIFIED AS A PROTEIN TARGET."""
<div align=Center>
  <a href="https://github.com/walla-42/Gene_Search">
    <img src="Images/Terminal_1.png" alt="Example1" width="500" height="200">
  </a>
</div>

The program will install all necessary files and create all the necessary data to identify drug targets. The program also provides useful graphical data so that the user can understand more about the drug target and the machine learning model used to identify the most probable drug target. 

<div align=Center>
  <a href="https://github.com/walla-42/Gene_Search">
    <img src="Images/Terminal_2.png" alt="Logo" width="600" height="200">
  </a>
</div>

The program will create a machine learning algorithm to predict the viabilty of certain drugs. The ML algorithm can be accessed by the pickle file created by the program. 

<div align=center>
  <a href="https://github.com/walla-42/Gene_Search">
    <img src="Images/Terminal_3.png" alt="Logo" width="450" height="100">
  </a>
</div>

<!--_For more examples, please refer to the [Documentation](https://example.com)_-->
<p align="right">(<a href="#readme-top">back to top</a>)</p> 



<!-- ROADMAP -->
## Roadmap

- [ ] Finish adding all the original algorithms in the OOP format -- Working on this one now
- [ ] I want to be able to do a lot more analysis of proteins. I have added the ability to do some basic analysis of proteins in the OOP format program. 
- [ ] Add more statistical analysis based on data generated as well as data visualization
- [ ] Add ways to export data generated in a document for the user -- Working on this one now

See the [open issues](https://github.com/walla-42/Gene_Search/issues) for a full list of proposed features (and known issues).

<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- CONTRIBUTING -->
## Contributing

Contributions are what make the open source community such an amazing place to learn, inspire, and create. Any contributions you make are **greatly appreciated**.

If you have a suggestion that would make this better, please fork the repo and create a pull request. You can also simply open an issue with the tag "enhancement".
Don't forget to give the project a star! Thanks again!

1. Fork the Project
2. Create your Feature Branch (`git checkout -b feature/AmazingFeature`)
3. Commit your Changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to the Branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request

<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- LICENSE -->
## License

Distributed under the MIT License. See `LICENSE.txt` for more information.

<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- CONTACT -->
## Contact


Project Link: [https://github.com/walla-42/GeneAnalyzer](https://github.com/walla-42/GeneAnalyzer)

<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- ACKNOWLEDGMENTS -->
## Acknowledgments

* [Walla-42](https://github.com/walla-42)


<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- MARKDOWN LINKS & IMAGES -->
<!-- https://www.markdownguide.org/basic-syntax/#reference-style-links -->
[contributors-shield]: https://img.shields.io/github/contributors/walla-42/Gene_Search.svg?style=for-the-badge
[contributors-url]: https://github.com/walla-42/Gene_Search/graphs/contributors
[forks-shield]: https://img.shields.io/github/forks/walla-42/Gene_Search.svg?style=for-the-badge
[forks-url]: https://github.com/walla-42/Gene_Search/network/members
[stars-shield]: https://img.shields.io/github/stars/walla-42/Gene_Search.svg?style=for-the-badge
[stars-url]: https://github.com/walla-42/Gene_Search/stargazers
[issues-shield]: https://img.shields.io/github/issues/walla-42/Gene_Search.svg?style=for-the-badge
[issues-url]: https://github.com/walla-42/Gene_Search/issues
[license-shield]: https://img.shields.io/github/license/walla-42/Gene_Search.svg?style=for-the-badge
[license-url]: https://github.com/walla-42/Gene_Search/blob/master/LICENSE.txt
[linkedin-shield]: https://img.shields.io/badge/-LinkedIn-black.svg?style=for-the-badge&logo=linkedin&colorB=555
[linkedin-url]: https://linkedin.com/in/walla42
[product-screenshot1]: Images/Terminal_1.png
[product-screenshot2]: Images/Terminal_2.png
[product-screenshot3]: Images/Terminal_3.png
[Python.js]: https://img.shields.io/badge/python-3670A0?style=for-the-badge&logo=python&logoColor=ffdd54
[Python-url]: https://www.python.org/
[SQlite.js]: https://img.shields.io/badge/sqlite-%2307405e.svg?style=for-the-badge&logo=sqlite&logoColor=white
[SQlite-url]: https://www.sqlite.org
