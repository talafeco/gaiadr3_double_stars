#! /usr/bin/python3
# USAGE: eszlelonaplo_uploader.py <name of the csv file to upload>
# CSV format:
#local_image_file_1,local_image_file_2,local_image_file_3,local_image_file_4,objectName,observerName,observationResult,observationPlace,observerEmail,observationDate,instrumentType,instrumentDiameter,instrumentFocalLength,instrumentResultingFocalLength,instrumentFocalLengthExtension,instrumentMagnification,instrumentEyepiece,instrumentEyepieceFocalLength,instrumentFilter,description,note,seeing,transparency,temperature,photoCamType,photoCamISO,photoCamGain,photoCamExpo,photoCamFrames,comp11,comp12,comp1sep,comp1pa,comp11Brightness,comp12Brightness,component11color,component12color,comp21,comp22,comp2sep,comp2pa,comp21Brightness,comp22Brightness,component21color,component22color,comp32,comp32,comp3sep,comp3pa,comp31Brightness,comp32Brightness,component31color,component32color,comp42,comp42,comp4sep,comp4pa,comp41Brightness,comp42Brightness,component41color,component42color,fov,endOfObservation,adminSection
#ds-1.jpg,ds-2.jpg,ds-3.jpg,ds-4.jpg,STFA1,Talabér Gergely,POSITIVE,Budapest,talafeco@gmail.com,2026-01-06 21:00:00,Refractor,100,1000,1500,1.5,150,Barium,10,ND1,Beautiful view of the Great Red Spot.,Testing upload script,7,4,5,ASI 120MM,100,100,100,10,A,B,10,100,5.5,7.5,Kék,Sárga,A,C,10,100,5.5,7.5,Kék,Sárga,A,D,10,100,5.5,7.5,Kék,Sárga,A,E,10,100,5.5,7.5,Kék,Sárga,20,00:01:00,BINARY

## To be done:
# 1. Authentication - OK
# 2. Add more images to upload - OK
# 3. Extend with pair 2, 3 and 4 and color - OK


import csv
import requests
import json
import os
import time
import sys

# --- CONFIGURATION ---

# 1. Credentials
USERNAME = "talafeco@gmail.com"  # Your observerId
PASSWORD = "ubWGA#Nw1:"

# 2. File Settings
IMAGE_FOLDER = "."
CSV_FILE = sys.argv[1]

# 3. API Endpoints
LOGIN_URL = "https://eszlelonaplo.mcse.hu/astro-log/api/users/login"
UPLOAD_URL = "https://eszlelonaplo.mcse.hu/astro-log/api/images/authenticated/upload"
CREATE_URL = "https://eszlelonaplo.mcse.hu/astro-log/api/observations/authenticated/create"

def login_and_get_session():
    """
    Logs in and sets up the session WITHOUT forcing Content-Type globally.
    """
    session = requests.Session()

    # Global headers for ALL requests
    session.headers.update({
        "User-Agent": "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/143.0.0.0 Safari/537.36",
        "Referer": "https://eszlelonaplo.mcse.hu/login",
        "Origin": "https://eszlelonaplo.mcse.hu",
        "Accept": "application/json, text/plain, */*"
    })

    print(f"[AUTH] Attempting to log in as {USERNAME}...")

    payload = {
        "observerId": USERNAME,
        "password": PASSWORD
    }

    try:
        # Pass Content-Type ONLY for login
        login_headers = {"Content-Type": "application/json"}
        response = session.post(LOGIN_URL, json=payload, headers=login_headers)

        if response.status_code == 200:
            print("[AUTH] Login Successful!")
            return session
        else:
            print(f"[AUTH] Login Failed. Status: {response.status_code}")
            print(f"Response: {response.text}")
            return None
    except Exception as e:
        print(f"[AUTH] Exception during login: {e}")
        return None

def upload_image(session, filename):
    filepath = os.path.join(IMAGE_FOLDER, filename)
    if not os.path.exists(filepath):
        print(f"[ERROR] File not found: {filepath}")
        return None

    print(f"[INFO] Uploading image: {filename}...")

    metadata_payload = {"adminSection": "BINARY", "fieldName": "imgDrawPhoto"}

    multipart_data = {
        'file': (filename, open(filepath, 'rb'), 'image/jpeg'),
        'imageMetadata': ('blob', json.dumps(metadata_payload), 'application/json')
    }

    try:
        response = session.post(UPLOAD_URL, files=multipart_data)

        if response.status_code == 200:
            try:
                data = response.json()
                if isinstance(data, dict):
                    return data.get('fileName') or data.get('filename') or data.get('file')
                return data
            except:
                return response.text
        else:
            print(f"[ERROR] Upload failed. Status: {response.status_code}")
            return None
    except Exception as e:
        print(f"[ERROR] Upload Exception: {e}")
        return None

def process_csv(session):
    if not os.path.exists(CSV_FILE):
        print(f"[ERROR] CSV file not found: {CSV_FILE}")
        return

    print(f"[START] Reading CSV: {CSV_FILE}")

    with open(CSV_FILE, mode='r', encoding='utf-8-sig') as f:
        reader = csv.DictReader(f)
        success_count = 0
        row_number = 0

        for row in reader:
            row_number += 1
            print(f"\n--- Processing Row {row_number}: {row.get('objectName')} ---")

            # --- 1. Prepare Payload & Handle Images ---
            payload = row.copy()

            # Define mapping: (CSV Column Name) -> (API Field Name)
            image_map = [
                ('local_image_file_1', 'imgDrawPhoto'),
                ('local_image_file_2', 'imgDetailDraw1'),
                ('local_image_file_3', 'imgDetailDraw2'),
                ('local_image_file_4', 'imgDetailDraw3')
            ]

            upload_failed = False

            for local_col, api_field in image_map:
                # Remove the local column from payload so we don't send garbage to server
                local_file = payload.pop(local_col, None)

                # Logic: Only attempt upload if a file is actually listed in the CSV
                if local_file and local_file.strip():
                    server_name = upload_image(session, local_file)

                    if server_name:
                        payload[api_field] = server_name
                    else:
                        # If user provided a file but it failed, we SHOULD skip to avoid bad data
                        print(f"[WARN] Failed to upload '{local_file}' for field '{api_field}'. Skipping row.")
                        upload_failed = True
                        break
                else:
                    # If no file provided, just set field to None (or don't send it)
                    # payload[api_field] = None # Optional: Explicitly send null
                    pass

            if upload_failed:
                continue

            # --- 2. Handle Number Conversions ---

            # LIST A: Text fields (Ignore these completely)
            text_fields = [
                'objectName', 'observerName', 'observationPlace', 'observerEmail', 'observationDate',
                'instrumentType', 'instrumentEyepiece', 'instrumentFilter', 'description', 'note',
                'adminSection', 'comp11', 'comp12', 'comp21', 'comp22', 'comp31', 'comp32', 'comp41', 'comp42',
                'component11color', 'component12color', 'component21color', 'component22color',
                'component31color', 'component32color', 'component41color', 'component42color',
                'observationResult', 'endOfObservation', 'imgDrawPhoto', 'imgDetailDraw1', 'imgDetailDraw2', 'imgDetailDraw3'
            ]

            # LIST B: Float fields that need "," instead of "." (e.g. "10.5" -> "10,5")
            float_comma_fields = [
                "comp1sep", "comp1pa", "comp11Brightness", "comp12Brightness",
                "comp2sep", "comp2pa", "comp21Brightness", "comp22Brightness",
                "comp3sep", "comp3pa", "comp31Brightness", "comp32Brightness",
                "comp4sep", "comp4pa", "comp41Brightness", "comp42Brightness",
                'instrumentFocalLengthExtension'
            ]

            for key, value in payload.items():
                if key in text_fields:
                    continue

                if value:
                    # Handle comma fields
                    if key in float_comma_fields:
                        try:
                            payload[key] = float(str(value).replace(',', '.')) # Ensure it is a valid float for JSON?
                            # WAIT: If the API wants a STRING with a comma (e.g. "10,5"), use this:
                            # payload[key] = str(value).replace('.', ',')

                            # If the API wants a JSON NUMBER (e.g. 10.5), use this:
                            payload[key] = float(str(value).replace(',', '.'))
                        except:
                            pass
                    else:
                        # Handle Integers
                        try:
                            payload[key] = int(float(value))
                        except ValueError:
                            pass

            # --- 3. Create Observation ---
            try:
                response = session.post(CREATE_URL, json=payload)

                if response.status_code in [200, 201]:
                    print(f"[SUCCESS] Observation created! ID: {response.json().get('id')}")
                    success_count += 1
                else:
                    print(f"[ERROR] Create failed. Status: {response.status_code}")
                    print(f"Response: {response.text}")
            except Exception as e:
                print(f"[ERROR] Request Exception: {e}")

            time.sleep(1)

    print(f"\n[DONE] Finished. Successfully uploaded {success_count} observations.")

if __name__ == "__main__":
    my_session = login_and_get_session()

    if my_session:
        process_csv(my_session)
