from selenium import webdriver
from selenium.webdriver.support.ui import Select
import time
from datetime import datetime
#print time.gmtime()
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
#import org.openqa.selenium.JavascriptExecutor

ye,mo,da,ho,mi,se,q,w,e = time.gmtime() # current time

w_d = 19 # wish date
w_h = 14 ; w_m = 59 ; w_s = 59 # library reservation opens at 00:00 (UTC +9) which is corresponding to 15:00 at UTC +0, Becuase of the delay it starts at 1 second before
ID = '2017324020' # id
PWD = '930211' # password

print time.gmtime()

############ To reduce the stand by time, it's going to sleep ##########
if w_h-ho >= 1:
    sleep = (w_h-ho)*3600.
elif 0 <= w_h-ho < 1:
    sleep = (w_m-mi)*30.
else:
    sleep = (w_h-ho+24)*3600.

print sleep/3600., 'hours sleep'

time.sleep(sleep)

ye,mo,da,ho,mi,se,q,w,e = time.gmtime()

print da,'day', ho, 'h', mi,'min', se,'s', 'wake up'
#########################################################################


while True :
    year, mon, day, hour, min, sec, a, b, c = time.gmtime()
    #print day, hour, min
    if day <= w_d and hour <= w_h and min <= w_m and sec < w_s :
        if sec==10:
            print hour, min,'not yet'

    if day==w_d and hour==w_h and min==w_m and sec==w_s:

        star_t = time.time()

        browser = webdriver.Chrome('C:\usr\chromedriver_win32\chromedriver')
        browser.get('https://library.yonsei.ac.kr/relation/seat')

        id = browser.find_element_by_name('id')
        id.send_keys(ID)

        pwd = browser.find_element_by_name('password')
        pwd.send_keys(PWD)

        click = browser.find_element_by_xpath('//*[@id="login"]/fieldset/div[2]/p')  # login click
        click.click()

        reserve = browser.find_element_by_xpath('/html/body/div[1]/app-root/index/div[1]/div/ul/li[2]')  # click reservation
        reserve.click()

        element = browser.find_element_by_xpath(
            "/html/body/div[1]/app-root/index/div[2]/app-facility/div/table/tbody/tr[1]/td[1]/perfect-scrollbar/div/div[1]/div/div[8]") # select date
        browser.execute_script("arguments[0].click();", element)  # click hide option such as latest date

        reserve3 = browser.find_element_by_xpath(
            '/html/body/div[1]/app-root/index/div[2]/app-facility/div/table/tbody/tr[1]/td[2]/perfect-scrollbar/div/div[1]/div/div[1]')  # select place
        reserve3.click()

        browser.implicitly_wait(1)

        reserve4 = browser.find_element_by_xpath(
            '/html/body/div[1]/app-root/index/div[2]/app-facility/div/table/tbody/tr[1]/td[3]/perfect-scrollbar/div/div[1]/div/div[3]')  # select floor
        reserve4.click()

        reserve5 = browser.find_element_by_xpath(
            '/html/body/div[1]/app-root/index/div[2]/app-facility/div/table/tbody/tr[1]/td[4]/perfect-scrollbar/div/div[1]/div/div[2]')  # select room
        reserve5.click()

        reserve6 = browser.find_element_by_xpath(
            '/html/body/div[1]/app-root/index/div[2]/app-facility/div/table/tbody/tr[1]/td[5]/perfect-scrollbar/div/div[1]/div/div[4]')  # select how long
        reserve6.click()

        reserve7 = browser.find_element_by_xpath(
            '/html/body/div[1]/app-root/index/div[2]/app-facility/div/table/tbody/tr[2]/td/div/div[2]/div[2]/div[1]/div[2]/div[3]')  # time select / last number 3 = 10am
        reserve7.click()

        browser.find_element_by_id("mat-select-4").click(); # click the purpose
        browser.find_element_by_css_selector('#mat-option-19').click()

        reason = browser.find_element_by_xpath(
            '/html/body/div[1]/app-root/index/div[2]/app-facility/div/table/tbody/tr[2]/td/div/div[2]/div[2]/div[3]/div/form/div[1]/div[2]/div[2]/div[1]/div[2]/div[2]/input')
        reason.send_keys('group meeting') # fill in the purpose

        id = browser.find_element_by_name('id')
        id.send_keys('2015313059')
        name = browser.find_element_by_name('name')
        name.send_keys('8514')
        click = browser.find_element_by_xpath(
            '/html/body/div[1]/app-root/index/div[2]/app-facility/div/table/tbody/tr[2]/td/div/div[2]/div[2]/div[3]/div/form/div[1]/div[2]/div[2]/div[2]/div/div/button')
        click.click()
        time.sleep(5)

        id = browser.find_element_by_name('id')
        id.send_keys('2012314108')
        name = browser.find_element_by_name('name')
        name.send_keys('2717')
        click = browser.find_element_by_xpath(
            '/html/body/div[1]/app-root/index/div[2]/app-facility/div/table/tbody/tr[2]/td/div/div[2]/div[2]/div[3]/div/form/div[1]/div[2]/div[2]/div[2]/div/div/button')
        click.click()
        time.sleep(5)

        id = browser.find_element_by_name('id')
        id.send_keys('2017323321')
        name = browser.find_element_by_name('name')
        name.send_keys('1030')
        click = browser.find_element_by_xpath(
            '/html/body/div[1]/app-root/index/div[2]/app-facility/div/table/tbody/tr[2]/td/div/div[2]/div[2]/div[3]/div/form/div[1]/div[2]/div[2]/div[2]/div/div/button')
        click.click()
        time.sleep(5)

        id = browser.find_element_by_name('id')
        id.send_keys('2017313062')
        name = browser.find_element_by_name('name')
        name.send_keys('5785')
        click = browser.find_element_by_xpath(
            '/html/body/div[1]/app-root/index/div[2]/app-facility/div/table/tbody/tr[2]/td/div/div[2]/div[2]/div[3]/div/form/div[1]/div[2]/div[2]/div[2]/div/div/button')
        click.click()
        time.sleep(5)

        submit = browser.find_element_by_xpath(
            "/html/body/div[1]/app-root/index/div[2]/app-facility/div/table/tbody/tr[2]/td/div/div[2]/div[2]/div[3]/div/form/div[2]/button")
        submit.click()

        f_t= time.time()

        print f_t-star_t
        print time.gmtime()

        print('reservation clear')
        break






