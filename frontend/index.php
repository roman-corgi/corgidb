<?php 
include "templates/header.php"; 
include "templates/headerclose.php"; 
?>

<h2> General Query </h2>
<p>
<form action="index.php" method="POST">
<textarea name="querytext" rows="12" cols="100">
<?php 
if (!empty($_GET["querytext"])){
    $sql =$_GET["querytext"];}
elseif (!empty($_POST["querytext"])){
    $sql = $_POST["querytext"]; }
else {$sql = 'SELECT main_id, spectype from Stars where sy_caltype = "RefStar"';}
    
echo "$sql";
?>

</textarea><br>
<!--<input type="submit">--!>

<button type="submit" name="submit" formaction="index.php">Submit</button>
<!-- <button type="submit" name="save" formaction="plansave.php">Save to CSV</button> --!>

</form>
</p>

<?php include("config.php"); ?>

<?php
$conn = new mysqli($servername, $username, $password, $dbname);
// Check connection
if ($conn->connect_error) {
    die("Connection failed: " . $conn->connect_error);
} 

echo "<h4>Query</h4><p>".$sql.";</p>\n\n";
$result = $conn->query($sql);
if ($result){
    if ($result->num_rows > 0) {
        echo "<h4>Result</h4><p>".$result->num_rows." rows returned.</p>\n\n";
        while($row = $result->fetch_assoc()) {
           if (empty($counter)){
              $counter = 0;
              echo "<div class='results-outer'>\n";
              echo "<table class='results' id='gentable'><thead><tr>\n";
              foreach(array_keys($row) as $paramName)
                echo "<th>".$paramName."</th>";
              echo "</tr></thead>\n";
           } 
           echo "<tr>";
           foreach(array_keys($row) as $paramName) {
               echo "<td>";
               if ($paramName == 'pl_name'){
                    echo "<a href='plandetail.php?name=".urlencode($row[$paramName]).
                    //"&scenario=".urlencode($row['scenario_name']).
                    //"&of_id=".urlencode($row['orbitfit_id']).
                    "'>".$row[$paramName]."</a>";                
               } else{
                   if (is_numeric($row[$paramName])){
                       echo number_format((float)$row[$paramName], 2, '.', '');

                   } else{
                       echo $row[$paramName];
                   }
               }
               echo "</td>";    
           }
           echo "</tr>";
         $counter++;
        }
        echo "</table></div>\n";
    }
    $result->close();
} else{
    echo "Query Error:\n".$conn->error;
}
$conn->close();
?>


<?php include "templates/footer.php"; ?>

