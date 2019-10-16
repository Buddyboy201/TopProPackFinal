a = '''import mechanize

br = mechanize.Browser()

br.set_handle_robots(False)

br.addheaders = [('User-agent', 'Firefox')]

br.open("http://dunbrack.fccc.edu/Guoli/PISCES_InputB.php")
print(br.geturl())

for f in br.forms():
    print(f)
print(br.geturl())
br.select_form(nr=0)
#print(br.form)
#br.form.add_file(open("C:\\alpha\\pdb_list.txt"), "text/plain", name="upload_file")
br["pasted"] = "1a91\n6qm4"
#print(br.form)
req = br.submit()

def select_form(form):
  return form.attrs.get('action', None) == 'PISCES_TakeUserInfo.php'

print("\n----------------\n")
#print(req.geturl())
#print(br.geturl())
#for f in br.forms():
#    print(f)
#br.form = br.global_form()
#br.select_form(action="http://dunbrack.fccc.edu/Guoli/PISCES_TakeUserInfo.php")
br.select_form(predicate=select_form)
print(br.form)
#br["userid"] = "Anish Prakriya"
#br["email"] = "researchprojecttest2052@gmail.com"
#br["institution"] = "UC Davis"
req2 = br.submit()
print(req2.geturl())'''

from selenium import webdriver

webpage = r"http://dunbrack.fccc.edu/Guoli/PISCES_InputB.php"

driver = webdriver.Chrome()
driver.get(webpage)
upload_file = driver.find_element_by_name("upload_file")
upload_file.send_keys("C:\\alpha\\pdb_list.txt")
submit = driver.find_element_by_name(".submit")
submit.click()
