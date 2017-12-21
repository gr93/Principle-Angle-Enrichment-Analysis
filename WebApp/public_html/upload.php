<?php
$target_dir = "uploads/";
$target_file1 = $target_dir . basename($_FILES["fileToUpload1"]["name"]);
$target_file2 = $target_dir . basename($_FILES["fileToUpload2"]["name"]);
$target_file3 = $target_dir . basename($_FILES["fileToUpload3"]["name"]);

$uploadOk = 1;
$imageFileType1 = pathinfo($target_file1,PATHINFO_EXTENSION);
$imageFileType2 = pathinfo($target_file2,PATHINFO_EXTENSION);
$imageFileType3 = pathinfo($target_file3,PATHINFO_EXTENSION);
$successcount = 0;

// Check if $uploadOk is set to 0 by an error
if ($uploadOk == 0) {
    echo "Sorry, your file was not uploaded.";
// if everything is ok, try to upload file
} else {
    if (move_uploaded_file($_FILES["fileToUpload1"]["tmp_name"], $target_file1)) {
        echo "The file ". basename( $_FILES["fileToUpload1"]["name"]). " has been uploaded.";
        $successcount = $successcount + 1;
    } else {
        echo "Sorry, there was an error uploading your expression data.";
    }
}
echo "<br />\n";
if ($uploadOk == 0) {
    echo "Sorry, your file was not uploaded.";
// if everything is ok, try to upload file
} else {
    if (move_uploaded_file($_FILES["fileToUpload2"]["tmp_name"], $target_file2)) {
        echo "The file ". basename( $_FILES["fileToUpload2"]["name"]). " has been uploaded.";
        $successcount = $successcount + 1;
    } else {
        echo "Sorry, there was an error uploading your sample classes.";
    }
}
echo "<br />\n";
if ($uploadOk == 0) {
    echo "Sorry, your file was not uploaded.";
// if everything is ok, try to upload file
} else {
    if (move_uploaded_file($_FILES["fileToUpload3"]["tmp_name"], $target_file3)) {
        echo "The file ". basename( $_FILES["fileToUpload3"]["name"]). " has been uploaded.";
        $successcount = $successcount + 1;
        if ($successcount = 3){
	    header('Location: https://web.njit.edu/~gr93/Rplots.pdf');
	}   
    } else {
        echo "Sorry, there was an error uploading your gamma values.";
    }
}
echo "<br />\n";
echo "<br />\n";
exec("cp uploads/expression_data.RData /afs/cad.njit.edu/u/g/r/gr93/R/x86_64-unknown-linux-gnu-library/3.1/GeoDE/data/ 2>&1");
exec("cp uploads/sampleclass.RData /afs/cad.njit.edu/u/g/r/gr93/R/x86_64-unknown-linux-gnu-library/3.1/GeoDE/data/");
exec("cp uploads/gammas.RData /afs/cad.njit.edu/u/g/r/gr93/R/x86_64-unknown-linux-gnu-library/3.1/GeoDE/data/");
echo "<br />\n";
exec("sh PAEA.sh 2>&1", $output);
foreach ($output as $line) {
    echo "$line\n";
    }
echo "<br />\n";
#header('Location: analysis.html');
?>
