<?php 
header("Content-Type:application/json");
include("config.php"); ?>
<?php

// Check inputs 
if (!isset($_GET['st_name'])){
    die("st_name not set");
}

$conn = new mysqli($servername, $username, $password, $dbname);
// Check connection
if ($conn->connect_error) {
    die("Connection failed: " . $conn->connect_error);
}

$st_name = $_GET['st_name'];
$sql = "SELECT `main_id` from StarAliases where alias = '$st_name'";

$result = $conn->query($sql);
$data = $result->fetch_all();

// Return JSON
echo json_encode($data);

$conn->close();
?>

