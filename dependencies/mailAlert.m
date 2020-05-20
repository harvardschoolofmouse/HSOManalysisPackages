function mailAlert(alertName, msg, sendto)
    % 
    %   Github HSOM version
    % 
    %   Created     5/1/2020 Allison Hamilos     ahamilos{at}g.harvard.edu
    %   Modified    5/1/2020 Allison Hamilos     ahamilos{at}g.harvard.edu
    % 
    if nargin < 3
        sendto = 'INSERT_YOUR_EMAIL_ADDRESSHERE'; %e.g., harvardschoolofmouse@HSOM.com
    end
    if strcmp(sendto, 'INSERT_YOUR_EMAIL_ADDRESSHERE')
        warning('mailAlert not configured. You can configure this code to mail you alerts when key analyses are complete by opening mailAlert.m')
    
    else
        if nargin < 2
            msg = [];
        end
        if nargin < 1
            alertName = 'none';
        end
        
            
        setpref('Internet','SMTP_Server','smtp.gmail.com'); % if your client is not gmail, you have to change this
        setpref('Internet','E_mail','########INSERTYOURHOSTEMAILADDRESS#########'); %e.g., HSOMhost@gmail.com
        setpref('Internet','SMTP_Username','########INSERTYOURHOSTaccount#########'); %e.g., HSOMhost
        setpref('Internet','SMTP_Password','########INSERTYOURPASSWORD#########'); %e.g., miceRgreat
        props = java.lang.System.getProperties;
        props.setProperty('mail.smtp.auth','true');
        props.setProperty('mail.smtp.socketFactory.class', 'javax.net.ssl.SSLSocketFactory');
        props.setProperty('mail.smtp.socketFactory.port','465');

        subject = ['Host:' getenv('computername') ' | ' alertName ' | ' datestr(now)];
        text = sprintf(['Host:' getenv('computername') ' | ' alertName ' | ' datestr(now) '\n\n' msg]);
        sendmail(sendto,subject, text);
        disp(['Alert: ' alertName ' mailed to ' sendto])
    end
end