<?php 
header("Content-Type:application/json");
include("config.php"); ?>
<?php
$conn = new mysqli($servername, $username, $password, $dbname);
// Check connection
if ($conn->connect_error) {
    die("Connection failed: " . $conn->connect_error);
} 

$sql = 'SELECT main_id, spectype, sy_vmag, sy_imag from Stars where sy_caltype = "RefStar"';
$result = $conn->query($sql);
$data = $result->fetch_all();

// Return JSON
echo json_encode($data);

$conn->close();
?>

