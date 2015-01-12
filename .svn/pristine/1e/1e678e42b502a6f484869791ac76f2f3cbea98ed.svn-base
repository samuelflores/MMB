#!/usr/bin/perl -w

# (C) Daniel Larsson 2012

@id = ('A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z','0','1','2','3','4','5','6','7','8','9','a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z');

$nr = -1;
$prev_id = "";
$previous_res = -1;
$first = 1;

while($line = <>) {

   if( $line =~ /^TER/ ){
       $nr++;

       if( $nr >= $#id ){
           $nr = 0;
       }
       printf $line;
       $previous_res = -1;
   }
   elsif( $line =~ m/(ATOM.{17})(.)(.{4})(.*)/ ){

       $pre = $1;
       $current_id = $2;
       $current_res = $3;
       $post = $4;

       $current_res =~ s/\s*(\d+)/$1/;

       if( $prev_id ne $current_id || $current_res < $previous_res || $current_res > $previous_res + 1 ) {
           $nr++;

           if( $nr >= $#id ){
               $nr = 0;
           }

           if( $first != 1 ){
               printf "TER\n";
           }
       }

       printf "%s%1s%4d%s\n", $pre, $id[$nr], $current_res, $post;

       if( $first == 1 ){
           $first = 0;
       }

       $prev_id = $current_id;
       $previous_res = $current_res;
   }
   elsif ( $line !~ /^TER/ ) {
       print $line;
   }

}
