# Aryl hydrocarbon receptor
The web page application for activity prediction for AhR (aryl hydrocarbon receptor)

![Logo AhR](https://github.com/nczub/AhR_Streamlit/blob/main/AhR_logo.png)

## Online Access
You can find the application at the following link:
[https://arylhydrocarbon-receptor.streamlit.app/](https://arylhydrocarbon-receptor.streamlit.app/)

## Docker

## Installation for Local calculation

To run the application locally, follow these steps:
1. **Clone the repository**
2. **Install dependencies**:

The needed packages are in file enviroment.txt

During installation you create conda environment 'for_AhR'

3. **Activate environment**
   
In the console activate conda environment:

```bash
$ conda activate for_AhR
```

4. Now, you can run the application:

```bash   
$ streamlit run Home.py
```
App should open in the browser or it will be available at 'http://localhost:8501'.

5. Finally, have fun and test my app!


## Batch mode

In batch mode, you can calculate predictions for multiple molecules. The online version of AhR has limitations based on Streamlit Cloud. The local app is much better to use for larger database.

Please, remember to upload CSV file with the column names 'smiles', because based on this system will predict class of activity.

## Author

I'm Natalia Łapińska and I'm the author of Aryl hydrocarbon receptor app and also SerotoninAI app [https://serotoninai.streamlit.app/](https://serotoninai.streamlit.app/)

- GitHub: https://github.com/nczub

- LinkedIn: https://www.linkedin.com/in/natalia-czub/

I would love to work with you if you want to collaborate on creating new QSAR or QSPR models.

## License

This project is available under the GNU General Public License v3.0 (GPL-3.0).
