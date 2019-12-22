#include <ttkCinemaQuery.h>

#include <vtkDataObject.h> // For port info
#include <vtkObjectFactory.h> // for new macro

#include <vtkDelimitedTextReader.h>
#include <vtkFieldData.h>
#include <vtkInformationVector.h>
#include <vtkSmartPointer.h>
#include <vtkTable.h>

#include <ttkUtils.h>

#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/replace.hpp>

vtkStandardNewMacro(ttkCinemaQuery);

ttkCinemaQuery::ttkCinemaQuery() {
  this->SetNumberOfInputPorts(5);
  this->SetNumberOfOutputPorts(1);
}
ttkCinemaQuery::~ttkCinemaQuery() {
}

int ttkCinemaQuery::FillInputPortInformation(int port, vtkInformation *info) {
  if(port < 0 || port > 4)
    return 0;

  info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkTable");
  info->Set(vtkAlgorithm::INPUT_IS_OPTIONAL(), 1);

  return 1;
}

int ttkCinemaQuery::FillOutputPortInformation(int port, vtkInformation *info) {
  if(port == 0)
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkTable");
  else
    return 0;
  return 1;
}

int ttkCinemaQuery::RequestData(vtkInformation *request,
                                vtkInformationVector **inputVector,
                                vtkInformationVector *outputVector) {
  ttk::Timer timer;

  // get output table
  auto outTable = vtkTable::GetData(outputVector);

  // ===========================================================================
  // Convert Input Tables to SQL Tables
  auto firstTable = vtkTable::GetData(inputVector[0]);

  std::vector<std::string> sqlTableDefinitions;
  std::vector<std::string> sqlInsertStatements;
  {
    ttk::Timer conversionTimer;
    this->printMsg("Converting input VTK tables to SQL tables", 0);

    for(size_t port = 0; port < 5; port++) {
      auto inVector = inputVector[port];
      if(inVector->GetNumberOfInformationObjects() != 1)
        continue;

      auto inTable = vtkTable::GetData(inputVector[port]);

      size_t nc = inTable->GetNumberOfColumns();
      size_t nr = inTable->GetNumberOfRows();
      std::vector<bool> isNumeric(nc);

      // -----------------------------------------------------------------------
      // Table Definition
      std::string sqlTableDefinition
        = "CREATE TABLE vtkTable" + std::to_string(port) + " (";
      for(size_t i = 0; i < nc; i++) {
        auto c = inTable->GetColumn(i);
        isNumeric[i] = c->IsNumeric();
        sqlTableDefinition += (i > 0 ? "," : "") + std::string(c->GetName())
                              + " " + (isNumeric[i] ? "REAL" : "TEXT");
      }
      sqlTableDefinition += ")";
      sqlTableDefinitions.push_back(sqlTableDefinition);

      // -----------------------------------------------------------------------
      // Insert Statements
      size_t q = 0;
      while(q < nr) {
        std::string sqlInsertStatement
          = "INSERT INTO vtkTable" + std::to_string(port) + " VALUES ";
        for(size_t j = 0; j < 500 && q < nr; j++) {
          if(j > 0)
            sqlInsertStatement += ",";

          sqlInsertStatement += "(";
          for(size_t i = 0; i < nc; i++) {
            if(i > 0)
              sqlInsertStatement += ",";
            if(isNumeric[i])
              sqlInsertStatement += inTable->GetValue(q, i).ToString();
            else
              sqlInsertStatement
                += "'" + inTable->GetValue(q, i).ToString() + "'";
          }
          sqlInsertStatement += ")";
          q++;
        }
        sqlInsertStatements.push_back(sqlInsertStatement);
      }
    }

    this->printMsg("Converting input VTK tables to SQL tables", 1,
                   conversionTimer.getElapsedTime(),
                   ttk::debug::LineMode::REPLACE);
  }

  // ===========================================================================
  // Replace Variables in Querystd::string (e.g. {time[2]})
  std::string finalQueryString;
  {
    std::string errorMsg;
    if(!ttkUtils::replaceVariables(this->GetSQLStatement(),
                                   firstTable->GetFieldData(), finalQueryString,
                                   errorMsg)) {
      this->printErr(errorMsg);
      return 0;
    }
  }

  // ===========================================================================
  // Compute Query Result
  std::stringstream csvResult;
  int csvNColumns = 0;
  int csvNRows = 0;

  int status
    = this->execute(sqlTableDefinitions, sqlInsertStatements, finalQueryString,
                    csvResult, csvNColumns, csvNRows);

  // ===========================================================================
  // Process Result
  {
#if VTK_MAJOR_VERSION < 7
    this->printMsg("This filter requires at least VTK 7.0", Priority::Error);
    return 0;
#else
    if(csvNRows < 1 || status != 1) {
      csvResult << "NULL";
      for(int i = 0; i < csvNColumns; i++)
        csvResult << ",NULL";
    }

    ttk::Timer conversionTimer;

    this->printMsg("Converting SQL result to VTK table", 0);

    auto reader = vtkSmartPointer<vtkDelimitedTextReader>::New();
    reader->SetReadFromInputString(true);
    reader->SetInputString(csvResult.str().data());
    reader->DetectNumericColumnsOn();
    reader->SetHaveHeaders(true);
    reader->SetFieldDelimiterCharacters(",");
    reader->Update();

    outTable->ShallowCopy(reader->GetOutput());

    auto outFD = outTable->GetFieldData();
    for(size_t port = 0; port < 5; port++) {
      auto inTable = vtkTable::GetData(inputVector[port]);
      if(!inTable)
        continue;
      auto inFD = inTable->GetFieldData();

      size_t n = inFD->GetNumberOfArrays();
      for(size_t i = 0; i < n; i++) {
        auto iArray = inFD->GetAbstractArray(i);
        if(!outFD->GetAbstractArray(iArray->GetName())) {
          outFD->AddArray(iArray);
        }
      }
    }

    this->printMsg("Converting SQL result to VTK table", 1,
                   conversionTimer.getElapsedTime(),
                   ttk::debug::LineMode::REPLACE);

#endif
  }

  // print stats
  this->printMsg(ttk::debug::Separator::L2);
  this->printMsg("Complete (#rows: " + std::to_string(csvNRows) + ")", 1,
                 timer.getElapsedTime());
  this->printMsg(ttk::debug::Separator::L1);

  return 1;
}
