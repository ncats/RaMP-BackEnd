from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import List


@dataclass
class SourceData:
    sourceId: str
    rampId: str
    IDtype: str
    geneOrCompound: str
    commonName: str
    priorityHMDBStatus: str
    dataSource: str

    def get_insert_format(self):
        return f"{self.sourceId}\t{self.rampId}\t{self.IDtype}\t{self.geneOrCompound}\t{self.commonName}\t{self.priorityHMDBStatus}\t{self.dataSource}"


class Analyte(ABC):
    rampId: str = ""
    representativeName: str = ""

    @abstractmethod
    def get_type(self):
        raise NotImplementedError('derived classes must implement this method')

    def get_insert_format(self):
        return f"{self.rampId}\t{self.get_type()}\t{self.representativeName}"

    @abstractmethod
    def getSourceData(self) -> List[SourceData]:
        raise NotImplementedError('derived classes must implement this method')

    def get_most_common_value(self, field="commonName", prefix=None):
        name_counts = {}
        source_info = self.getSourceData()
        for source_data in source_info:
            value = getattr(source_data, field)
            if value is None or value == 'None' or value == 'NA':
                continue
            if prefix is not None and not value.startswith(prefix):
                continue
            lowercase_value = value.lower()
            if lowercase_value in name_counts:
                name_counts[lowercase_value]['count'] += 1
            else:
                name_counts[lowercase_value] = {
                    'value': value,
                    'count': 1
                }
        sorted_names = dict(sorted(name_counts.items(), key=lambda item: item[1]['count'], reverse=True))

        best_names = []
        highest_count = None
        for lowercase_value, name_data in sorted_names.items():
            count = name_data['count']
            value = name_data['value']
            if highest_count is None:
                highest_count = count
            if count == highest_count:
                best_names.append(value)

        if len(best_names) == 0:
            if prefix is not None:
                return None
            if field == 'commonName':
                return self.get_most_common_value('sourceId')
            else:
                return self.rampId

        return best_names[0]
