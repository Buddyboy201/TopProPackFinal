
import imaplib
import email

# Your IMAP Settings
host = 'imap.gmail.com'
user = 'researchprojecttest2052@gmail.com'
password = 'cmVzZWFyY2hwcm9qZWN0dGVzdA=='


def getLink(msg, start, end):
    return (str(msg)[start:end]).strip()


# Connect to the server
print('Connecting to ' + host)
mailBox = imaplib.IMAP4_SSL(host)

# Login to our account
mailBox.login(user, password)

mailBox.select()
searchQuery = '(FROM "Roland.Dunbrack@fccc.edu")'

result, data = mailBox.uid('search', None, searchQuery)
ids = data[0]
# list of uids
id_list = ids.split()

i = len(id_list)
#for x in range(1):
latest_email_uid = id_list[i-1]
result, email_data = mailBox.uid('fetch', latest_email_uid, '(RFC822)')
raw_email = email_data[0][1]
raw_email_string = raw_email.decode('utf-8')
email_message = email.message_from_string(raw_email_string)

s1 = str(email_message).find('And check your original input file:')
s2 = str(email_message).find('Download sequence ID list file')
s3 = str(email_message).find('Download fasta format sequence file')
s4 = str(email_message).find('Download similarity log file')
s5 = str(email_message).find('Download all detected pairwise scores')
s6 = str(email_message).find('________________________________')
input_file = getLink(email_message, s1+len('And check your original input file:'), s2)
sequence_id_file = getLink(email_message, s2+len('Download sequence ID list file'), s3)
fasta_file = getLink(email_message, s3+len('Download fasta format sequence file'), s4)
similarity_log_file = getLink(email_message, s4+len('Download similarity log file'), s5)
pairwise_scores_file = getLink(email_message, s5+len('Download all detected pairwise scores'), s6)
print(input_file)
print(sequence_id_file)
print(fasta_file)
print(similarity_log_file)
print(pairwise_scores_file)
mailBox.close()
mailBox.logout()
