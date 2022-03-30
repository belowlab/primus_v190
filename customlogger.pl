use Term::ANSIColor qw(:constants);
use strict;
use warnings;

# -------------------------------------
# USAGE : require './customlogger.pl';
# -------------------------------------

sub custom_log {

    my $message = shift;
    my $message_type = shift;
    my $color;

    my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
    my $timestamp = sprintf ( "%02d:%02d:%02d",$hour,$min,$sec);
    
    if ($message_type eq 'success') {
        print GREEN, "[$timestamp] $message\n", RESET;
    }
    if ($message_type eq 'error') {
        print RED, "[$timestamp] $message\n", RESET;
    }
    if ($message_type eq 'debug') {
        print YELLOW, "[$timestamp] $message\n", RESET;
    }
    if ($message_type eq 'info') {
        print BLUE, "[$timestamp] $message\n", RESET;
    }

}

return 1;