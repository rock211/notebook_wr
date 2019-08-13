from selenium import webdriver
from selenium.webdriver.support.ui import Select
import time
from datetime import datetime
#print time.gmtime()

ye,mo,da,ho,mi,se,q,w,e = time.gmtime() # current time

ID = '2017324020' # id
PWD = '930211' # password
DATE = '20190311' # choose the wish date
w_d = 9 # wish date - 2 days = reservation date
w_h = 4 ; w_m = 59 ; w_s = 59 # shuttle bus reservation opens at 2 pm (UTC +9) which is corresponding to 5 am in UTC +0, Becuase of the delay it starts at 1 second before

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
###########################################################################

while True :
    year, mon, day, hour, min, sec, a, b, c = time.gmtime()
    #print day, hour, min
    if day <= w_d and hour <= w_h and min <= w_m and sec < w_s :
        if sec==10:
            print hour, min,'not yet'

    if day==w_d and hour==w_h and min==w_m and sec==w_s:

        star_t = time.time()

        ############ from inter to shinchon
        browser = webdriver.Chrome('C:\usr\chromedriver_win32\chromedriver')
        browser.get('https://ysweb.yonsei.ac.kr/ysbus.jsp') #website

        id = browser.find_element_by_name('userid') #id
        id.send_keys(ID)

        pwd = browser.find_element_by_name('password') #password
        pwd.send_keys(PWD)

        click = browser.find_element_by_class_name('submit') #login
        click.click()

        reserve = browser.find_element_by_partial_link_text("(RESERVATION)") # click reservation
        reserve.click()

        date = Select(browser.find_element_by_id("selectedDate")) # select date
        date.select_by_value(DATE)

        posi = browser.find_element_by_xpath('//div/ul[2]/li[2]/a/span')  # choose international campus
        posi.click()

        cell = browser.find_element_by_xpath('//table/tbody/tr[9]/td[5]/select/option[2]') # tr[n] means the order of the wish time
        cell.click()
        confirm = browser.find_element_by_xpath('//table/tbody/tr[9]/td[6]/a') # tr[n] means the order of the wish time
        confirm.click()

        ########## from shinchon to inter

        browser = webdriver.Chrome('C:\usr\chromedriver_win32\chromedriver')
        browser.get('https://ysweb.yonsei.ac.kr/ysbus.jsp')

        id = browser.find_element_by_name('userid')
        id.send_keys(ID)

        pwd = browser.find_element_by_name('password')
        pwd.send_keys(PWD)

        click = browser.find_element_by_class_name('submit')
        click.click()

        reserve = browser.find_element_by_partial_link_text("(RESERVATION)")
        reserve.click()

        date = Select(browser.find_element_by_id("selectedDate"))
        date.select_by_value(DATE) # select the data

        cell = browser.find_element_by_xpath('//table/tbody/tr[5]/td[5]/select/option[2]')# tr[n] means the order of the wish time
        cell.click()

        confirm = browser.find_element_by_xpath('//table/tbody/tr[5]/td[6]/a') # tr[n] means the order of the wish time
        confirm.click()

        f_t= time.time()

        print f_t-star_t
        print time.gmtime()

        print('reservation clear')
        break



'''
browser = webdriver.Chrome('C:\usr\chromedriver_win32\chromedriver')
browser.get('https://ysweb.yonsei.ac.kr/ysbus.jsp')

id = browser.find_element_by_name('userid')
id.send_keys('2017324020')

pwd = browser.find_element_by_name('password')
pwd.send_keys('930211')

click = browser.find_element_by_class_name('submit')
click.click()

reserve = browser.find_element_by_partial_link_text("(RESERVATION)")
reserve.click()

#desti = browser.find_element_by_class_name("on")
#desti.click()

date = Select(browser.find_element_by_id("selectedDate"))
date.select_by_value("20190308") # select the data

#posi = browser.find_element_by_xpath('//div/ul[2]/li[2]/a/span') # choose international campus
#posi.click()

# tr : 12 = 16:30 / 13 17:00

cell = browser.find_element_by_xpath('//table/tbody/tr[4]/td[5]/select/option[2]') # tr means the order of wish time
cell.click()

confirm = browser.find_element_by_xpath('//table/tbody/tr[4]/td[6]/a') # tr means the order of wish time
confirm.click()
'''
#year,mon,day,hour,min,sec, a,b,c = time.gmtime()
#print day, hour
