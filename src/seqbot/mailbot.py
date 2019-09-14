#!/usr/bin/env python

import pathlib
import smtplib
import subprocess

from email.message import EmailMessage


def config_msg(subject: str, sender: str, mailto: str, content: str):
    msg = EmailMessage()
    msg.set_content(content)
    msg["Subject"] = subject
    msg["From"] = sender
    msg["To"] = mailto

    return msg


def send_mail(msg: EmailMessage):
    with smtplib.SMTP("localhost") as smtp:
        smtp.send_message(msg)


def demux_mail(s3_uri: str, run_name: str, email_config: dict, no_mismatch: bool):
    if no_mismatch:
        header = "NOTE! This run was demuxed with --barcode_mismatches 0\n\n"
    else:
        header = ""

    msg = config_msg(
        subject=f"[Seqbot] demux for {run_name} is complete!",
        sender=email_config["username"],
        mailto=",".join(email_config["addresses_to_email"]),
        content=f"""{header}Results are located in:
    {s3_uri}/{run_name}

- seqbot
""",
    )

    send_mail(msg)


def mail_nova_index(run_name: str, index_counts: pathlib.Path, email_config: dict):
    msg = config_msg(
        subject=f"[Seqbot] index counts for {run_name}",
        sender=email_config["username"],
        mailto=",".join(email_config["addresses_to_email"]),
        content=f"""
The most common indexes are attached as a {index_counts.suffix} file.

- seqbot
""",
    )

    with open(index_counts, "rb") as f:
        if index_counts.suffix == ".gz":
            msg.add_attachment(
                f.read(), "text", "plain+gzip", filename=index_counts.name
            )
        else:
            msg.add_attachment(f.read(), "text", "plain", filename=index_counts.name)

    send_mail(msg)


def error_mail(run_name: str, proc: subprocess.CompletedProcess, email_config: dict):
    stderr_lines = proc.stderr.splitlines()

    n_lines = min(10, len(stderr_lines) // 2)

    head = "\n".join(stderr_lines[:n_lines])
    tail = "\n".join(stderr_lines[-n_lines:])

    msg = config_msg(
        subject=f"[Seqbot] demux for {run_name} had an error",
        sender=email_config["username"],
        mailto=",".join(email_config["addresses_to_email_on_error"]),
        content=f"""There was an error while demuxing run {run_name}:

{head}
...
{tail}

- seqbot
""",
    )

    msg.add_attachment(proc.stderr, filename=f"{run_name}_error.txt")

    send_mail(msg, email_config)


def samplesheet_error_mail(run_name: str, email_config: dict):
    msg = config_msg(
        subject=f"[Seqbot] {run_name} might have a bad samplesheet",
        sender=email_config["username"],
        mailto=",".join(email_config["addresses_to_email"]),
        content=f"""There were non-ASCII characters in the file, which is usually bad.

- seqbot
""",
    )
